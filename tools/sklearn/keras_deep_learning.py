import argparse
import json
import keras
import pandas as pd
import pickle
import six
import warnings

from ast import literal_eval
from keras.models import Sequential, Model
from galaxy_ml.utils import try_get_attr, get_search_params, SafeEval


safe_eval = SafeEval()


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
    if constraint_type in ('None', ''):
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
        if value in ('None', ''):
            params[key] = None
            continue

        if type(value) in [int, float, bool]\
                or (type(value) is str and value.isalpha()):
            continue

        if key in ['input_shape', 'noise_shape', 'shape', 'batch_shape',
                   'target_shape', 'dims', 'kernel_size', 'strides',
                   'dilation_rate', 'output_padding', 'cropping', 'size',
                   'padding', 'pool_size', 'axis', 'shared_axes'] \
                and isinstance(value, str):
            params[key] = _handle_shape(value)

        elif key.endswith('_regularizer') and isinstance(value, dict):
            params[key] = _handle_regularizer(value)

        elif key.endswith('_constraint') and isinstance(value, dict):
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
        kwargs = options.pop('kwargs', '')

        # parameters needs special care
        options = _handle_layer_parameters(options)

        if kwargs:
            kwargs = safe_eval('dict(' + kwargs + ')')
            options.update(kwargs)

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
        kwargs = options.pop('kwargs', '')

        # parameters needs special care
        options = _handle_layer_parameters(options)

        if kwargs:
            kwargs = safe_eval('dict(' + kwargs + ')')
            options.update(kwargs)

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
        json.dump(json.loads(json_string), f, indent=2)


def build_keras_model(inputs, outfile, model_json, infile_weights=None,
                      batch_mode=False, outfile_params=None):
    """ for `keras_model_builder` tool

    Parameters
    ----------
    inputs : dict
        loaded galaxy tool parameters from `keras_model_builder` tool.
    outfile : str
        Path to galaxy dataset containing the keras_galaxy model output.
    model_json : str
        Path to dataset containing keras model JSON.
    infile_weights : str or None
        If string, path to dataset containing model weights.
    batch_mode : bool, default=False
        Whether to build online batch classifier.
    outfile_params : str, default=None
        File path to search parameters output.
    """
    with open(model_json, 'r') as f:
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

        train_metrics = (inputs['mode_selection']['compile_params']
                         ['metrics']).split(',')
        if train_metrics[-1] == 'none':
            train_metrics = train_metrics[:-1]
        options['metrics'] = train_metrics

        options.update(inputs['mode_selection']['fit_params'])
        options['seed'] = inputs['mode_selection']['random_seed']

        if batch_mode:
            generator = get_batch_generator(inputs['mode_selection']
                                            ['generator_selection'])
            options['data_batch_generator'] = generator
            options['prediction_steps'] = \
                inputs['mode_selection']['prediction_steps']
            options['class_positive_factor'] = \
                inputs['mode_selection']['class_positive_factor']
        estimator = klass(config, **options)
        if outfile_params:
            hyper_params = get_search_params(estimator)
            # TODO: remove this after making `verbose` tunable
            for h_param in hyper_params:
                if h_param[1].endswith('verbose'):
                    h_param[0] = '@'
            df = pd.DataFrame(hyper_params, columns=['', 'Parameter', 'Value'])
            df.to_csv(outfile_params, sep='\t', index=False)

    print(repr(estimator))
    # save model by pickle
    with open(outfile, 'wb') as f:
        pickle.dump(estimator, f, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    warnings.simplefilter('ignore')

    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-m", "--model_json", dest="model_json")
    aparser.add_argument("-t", "--tool_id", dest="tool_id")
    aparser.add_argument("-w", "--infile_weights", dest="infile_weights")
    aparser.add_argument("-o", "--outfile", dest="outfile")
    aparser.add_argument("-p", "--outfile_params", dest="outfile_params")
    args = aparser.parse_args()

    input_json_path = args.inputs
    with open(input_json_path, 'r') as param_handler:
        inputs = json.load(param_handler)

    tool_id = args.tool_id
    outfile = args.outfile
    outfile_params = args.outfile_params
    model_json = args.model_json
    infile_weights = args.infile_weights

    # for keras_model_config tool
    if tool_id == 'keras_model_config':
        config_keras_model(inputs, outfile)

    # for keras_model_builder tool
    else:
        batch_mode = False
        if tool_id == 'keras_batch_models':
            batch_mode = True

        build_keras_model(inputs=inputs,
                          model_json=model_json,
                          infile_weights=infile_weights,
                          batch_mode=batch_mode,
                          outfile=outfile,
                          outfile_params=outfile_params)
