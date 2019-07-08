import keras
import six

from ast import literal_eval
from keras.models import Sequential, Model
from galaxy_ml.utils import try_get_attr


def _handle_shape(literal):
    """Eval integer or list/tuple of integers from string

    Parameters:
    -----------
    literal : str.
    """
    literal = literal.strip()
    if not literal:
        return None
    try:
        return literal_eval(literal)
    except NameError as e:
        print(e)
        return literal


def _handle_regularizer(literal):
    """Construct regularizer from string literal

    Parameters
    ----------
    literal : str. E.g. '(0.1, 0)'
    """
    literal = literal.strip()
    if not literal:
        return None

    l1, l2 = literal_eval(literal)

    if not l1 and not l2:
        return None

    if l1 is None:
        l1 = 0.
    if l2 is None:
        l2 = 0.

    return keras.regularizers.l1_l2(l1=l1, l2=l2)


def _handle_constraint(config):
    """Construct constraint from galaxy tool parameters.
    Suppose correct dictionary format

    Parameters
    ----------
    config : dict. E.g.
        "bias_constraint":
            {"constraint_options":
                {"max_value":1.0,
                "min_value":0.0,
                "axis":"[0, 1, 2]"
                },
            "constraint_type":
                "MinMaxNorm"
            }
    """
    constraint_type = config['constraint_type']
    if constraint_type == 'None':
        return None

    klass = getattr(keras.constraints, constraint_type)
    options = config.get('constraint_options', {})
    if 'axis' in options:
        options['axis'] = literal_eval(options['axis'])

    return klass(**options)


def _handle_lambda(literal):
    return None


def _handle_layer_parameters(params):
    """Access to handle all kinds of parameters
    """
    for key, value in six.iteritems(params):
        if value == 'None':
            params[key] = None
            continue

        if type(value) in [int, float, bool]\
                or (type(value) is str and value.isalpha()):
            continue

        if key in ['input_shape', 'noise_shape', 'shape', 'batch_shape',
                   'target_shape', 'dims', 'kernel_size', 'strides',
                   'dilation_rate', 'output_padding', 'cropping', 'size',
                   'padding', 'pool_size', 'axis', 'shared_axes']:
            params[key] = _handle_shape(value)

        elif key.endswith('_regularizer'):
            params[key] = _handle_regularizer(value)

        elif key.endswith('_constraint'):
            params[key] = _handle_constraint(value)

        elif key == 'function':  # No support for lambda/function eval
            params.pop(key)

    return params


def get_sequential_model(config):
    """Construct keras Sequential model from Galaxy tool parameters

    Parameters:
    -----------
    config : dictionary, galaxy tool parameters loaded by JSON
    """
    model = Sequential()
    input_shape = _handle_shape(config['input_shape'])
    layers = config['layers']
    for layer in layers:
        options = layer['layer_selection']
        layer_type = options.pop('layer_type')
        klass = getattr(keras.layers, layer_type)
        other_options = options.pop('layer_options', {})
        options.update(other_options)

        # parameters needs special care
        options = _handle_layer_parameters(options)

        # add input_shape to the first layer only
        if not getattr(model, '_layers') and input_shape is not None:
            options['input_shape'] = input_shape

        model.add(klass(**options))

    return model


def get_functional_model(config):
    """Construct keras functional model from Galaxy tool parameters

    Parameters
    -----------
    config : dictionary, galaxy tool parameters loaded by JSON
    """
    layers = config['layers']
    all_layers = []
    for layer in layers:
        options = layer['layer_selection']
        layer_type = options.pop('layer_type')
        klass = getattr(keras.layers, layer_type)
        inbound_nodes = options.pop('inbound_nodes', None)
        other_options = options.pop('layer_options', {})
        options.update(other_options)

        # parameters needs special care
        options = _handle_layer_parameters(options)
        # merge layers
        if 'merging_layers' in options:
            idxs = literal_eval(options.pop('merging_layers'))
            merging_layers = [all_layers[i-1] for i in idxs]
            new_layer = klass(**options)(merging_layers)
        # non-input layers
        elif inbound_nodes is not None:
            new_layer = klass(**options)(all_layers[inbound_nodes-1])
        # input layers
        else:
            new_layer = klass(**options)

        all_layers.append(new_layer)

    input_indexes = _handle_shape(config['input_layers'])
    input_layers = [all_layers[i-1] for i in input_indexes]

    output_indexes = _handle_shape(config['output_layers'])
    output_layers = [all_layers[i-1] for i in output_indexes]

    return Model(inputs=input_layers, outputs=output_layers)


def get_batch_generator(config):
    """Construct keras online data generator from Galaxy tool parameters

    Parameters
    -----------
    config : dictionary, galaxy tool parameters loaded by JSON
    """
    generator_type = config.pop('generator_type')
    klass = try_get_attr('galaxy_ml.preprocessors', generator_type)

    if generator_type == 'GenomicIntervalBatchGenerator':
        config['ref_genome_path'] = 'to_be_determined'
        config['intervals_path'] = 'to_be_determined'
        config['target_path'] = 'to_be_determined'
        config['features'] = 'to_be_determined'
    else:
        config['fasta_path'] = 'to_be_determined'

    return klass(**config)


def config_keras_model(inputs, outfile):
    """ config keras model layers and output JSON

    Parameters
    ----------
    inputs : dict
        loaded galaxy tool parameters from `keras_model_config`
        tool.
    outfile : str
        Path to galaxy dataset containing keras model JSON.
    """
    model_type = inputs['model_selection']['model_type']
    layers_config = inputs['model_selection']

    if model_type == 'sequential':
        model = get_sequential_model(layers_config)
    else:
        model = get_functional_model(layers_config)

    json_string = model.to_json()

    with open(outfile, 'w') as f:
        f.write(json_string)


def build_keras_model(inputs, outfile, infile_json, infile_weights=None,
                      batch_mode=False):
    """ for `keras_model_builder` tool

    Parameters
    ----------
    inputs : dict
        loaded galaxy tool parameters from `keras_model_builder`
        tool.
    outfile : str
        Path to galaxy dataset containing keras_galaxy model.
    infile_json : str
        Path to dataset containing keras model JSON
    infile_weights : str or None
        If string, path to dataset containing model weights.
    batch_mode : bool, default=False
        Whether to build online batch classifier.
    """
    with open(infile_json, 'r') as f:
        json_model = json.load(f)

    config = json_model['config']

    options = {}

    if json_model['class_name'] == 'Sequential':
        options['model_type'] = 'sequential'
        klass = Sequential
    elif json_model['class_name'] == 'Model':
        options['model_type'] = 'functional'
        klass = Model
    else:
        raise ValueError("Unknow Keras model class: %s"
                         % json_model['class_name'])

    # load prefitted model
    if inputs['mode_selection']['mode_type'] == 'prefitted':
        estimator = klass.from_config(config)
        estimator.load_weights(infile_weights)
    # build train model
    else:
        cls_name = inputs['mode_selection']['learning_type']
        klass = try_get_attr('galaxy_ml.keras_galaxy_models', cls_name)

        options['loss'] = (inputs['mode_selection']
                           ['compile_params']['loss'])
        options['optimizer'] =\
            (inputs['mode_selection']['compile_params']
             ['optimizer_selection']['optimizer_type']).lower()

        options.update((inputs['mode_selection']['compile_params']
                        ['optimizer_selection']['optimizer_options']))
        options.update(inputs['mode_selection']['fit_params'])

        if batch_mode:
            generator = get_batch_generator(inputs['mode_selection']
                                            ['generator_selection'])
            options['data_batch_generator'] = generator
        estimator = klass(config, **options)

    print(repr(estimator))
    # save model by pickle
    with open(outfile, 'wb') as f:
        pickle.dump(estimator, f, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    import json
    import pickle
    import sys
    import warnings

    warnings.simplefilter('ignore')
    input_json_path = sys.argv[1]
    with open(input_json_path, 'r') as param_handler:
        inputs = json.load(param_handler)

    tool_id = sys.argv[2]
    outfile = sys.argv[3]

    if len(sys.argv) > 4:
        infile_json = sys.argv[4]
    if len(sys.argv) > 5:
        infile_weights = sys.argv[5]
    else:
        infile_weights = None

    # for keras_model_builder tool
    if tool_id == 'keras_model_builder':
        build_keras_model(inputs, outfile, infile_json,
                          infile_weights=infile_weights)

    elif tool_id == 'keras_batch_models':
        build_keras_model(inputs, outfile, infile_json,
                          infile_weights=infile_weights,
                          batch_mode=True)

    # for keras_model_config tool
    else:
        config_keras_model(inputs, outfile)
