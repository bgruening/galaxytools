"""
Classes:

    ModelToDict
    DictToModel
    H5Model

Functions:

    dump(object) -> dictionary
    load(dictionary) -> object
    save_model(dictionary) -> h5
    load_model(h5) -> dictionary

"""

import sys
import six
import types
import numpy
import h5py
import json

# reserved keys
_PY_VERSION = '-cpython-'
_OBJ = '-object-'
_NP_VERSION = '-numpy-'
_REDUCE = '-reduce-'
_GLOBAL = '-global-'
_FUNC = '--func-'
_ARGS = '-args-'
_STATE = '-state-'
_KEYS = '-keys-'
_MEMO = '-memo-'
_NP_NDARRAY = '-np_ndarray-'
_DTYPE = '-dtype-'
_VALUES = '-values-'
_NP_DATATYPE = '-np_datatype-'
_DATATYPE = '-datatype-'
_VALUE = '-value-'
_TUPLE = '-tuple-'
_SET = '-set-'
_BYTES = '-bytes-'


class JsonPicklerError(Exception):
    pass


class ModelToDict:
    """
    Follow the track of python `pickle`
    Turn a scikit-learn model to a JSON-compatiable dictionary
    """
    def __init__(self):
        self.memo = {}

    def clear_memo(self):
        """
        Clears the `memo`
        """
        self.memo.clear()

    def memoize(self, obj):
        """
        Store an object id in the `memo`
        """
        assert id(obj) not in self.memo
        idx = len(self.memo)
        self.memo[id(obj)] = idx, obj

    def dump(self, obj):
        """
        Main access of object save
        """
        py_version = sys.version.split(' ')[0]
        retv = {_PY_VERSION: py_version}
        np_version = sys.modules.get('numpy').__version__
        if np_version:
            retv[_NP_VERSION] = np_version
        retv[_OBJ] = self.save(obj)
        return retv

    def save(self, obj):

        # Check the `memo``
        x = self.memo.get(id(obj))
        if x:
            rval = {_MEMO: x[0]}
            return rval

        # Check type in `dispath` table
        t = type(obj)
        f = self.dispatch.get(t)
        if f:
            return f(self, obj)

        # Check for a class with a custom metaclass; treat as regular class
        try:
            issc = issubclass(t, type)
        except TypeError:
            issc = 0
        if issc:
            return self.save_global(obj)

        return self.save_reduce(obj)

    def save_reduce(self, obj):
        """
        Decompose an object using pickle reduce
        """
        reduce = getattr(obj, "__reduce__", None)
        if reduce:
            rv = reduce()
        else:
            raise JsonPicklerError("Can't reduce %r object: %r" % (obj.__name__, obj))
        assert (type(rv) is tuple), \
            "%s must return a tuple, but got %s" % (reduce, type(rv))
        l = len(rv)
        assert (l in [2, 3]),\
            "Reduce tuple is expected to return 2- 3 elements, but got %d elements" % l

        save = self.save

        retv = {}

        func = rv[0]
        assert callable(func), "func from reduce is not callable"
        retv[_FUNC] = save(func)

        args = rv[1]
        assert (type(args) is tuple)
        retv[_ARGS] = {_TUPLE: save(list(args))}

        if l == 3:
            state = rv[2]
            retv[_STATE] = save(state)

        self.memoize(obj)
        return {_REDUCE: retv}

    dispatch = {}

    def save_primitive(self, obj):
        return obj

    dispatch[type(None)] = save_primitive
    dispatch[bool] = save_primitive
    dispatch[int] = save_primitive
    if six.PY2:
        dispatch[long] = save_primitive
    dispatch[float] = save_primitive
    dispatch[complex] = save_primitive

    def save_bytes(self, obj):
        self.memoize(obj)
        return {_BYTES: obj.decode('utf-8')}

    dispatch[bytes] = save_bytes

    def save_string(self, obj):
        self.memoize(obj)
        return obj

    dispatch[str] = save_string

    if six.PY2:
        dispatch[unicode] = save_string

    def save_list(self, obj):
        return [self.save(e) for e in obj]

    dispatch[list] = save_list

    def save_tuple(self, obj):
        aslist = self.save(list(obj))
        return {_TUPLE: aslist}

    dispatch[tuple] = save_tuple

    def save_set(self, obj):
        aslist = self.save(list(obj))
        return {_SET: aslist}

    dispatch[set] = save_set

    def save_dict(self, obj):
        newdict = {}
        _keys = list(obj.keys())
        _keys.sort()
        newdict[_KEYS] = _keys
        for k in _keys:
            newdict[k] = self.save(obj[k])
        return newdict

    dispatch[dict] = save_dict

    def save_global(self, obj):
        name = getattr(obj, '__name__', None)
        if name is None:
            raise JsonPicklerError("Can't get global name for object %r" % obj)
        module_name = getattr(obj, '__module__', None)
        if module_name is None:
            raise JsonPicklerError("Can't get global module name for object %r" % obj)

        newdict = {_GLOBAL: [module_name, name]}
        self.memoize(obj)
        return newdict

    dispatch[types.FunctionType] = save_global
    dispatch[types.BuiltinFunctionType] = save_global

    def save_np_ndarray(self, obj):
        newdict = {}
        newdict[_DTYPE] = self.save(obj.dtype)
        newdict[_VALUES] = self.save(obj.tolist())
        return {_NP_NDARRAY: newdict}

    dispatch[numpy.ndarray] = save_np_ndarray

    def save_np_datatype(self, obj):
        newdict = {}
        newdict[_DATATYPE] = self.save(type(obj))
        newdict[_VALUE] = self.save(obj.item())
        return {_NP_DATATYPE: newdict}

    dispatch[numpy.bool_] = save_np_datatype
    dispatch[numpy.int_] = save_np_datatype
    dispatch[numpy.intc] = save_np_datatype
    dispatch[numpy.intp] = save_np_datatype
    dispatch[numpy.int8] = save_np_datatype
    dispatch[numpy.int16] = save_np_datatype
    dispatch[numpy.int32] = save_np_datatype
    dispatch[numpy.int64] = save_np_datatype
    dispatch[numpy.uint8] = save_np_datatype
    dispatch[numpy.uint16] = save_np_datatype
    dispatch[numpy.uint32] = save_np_datatype
    dispatch[numpy.uint64] = save_np_datatype
    dispatch[numpy.float_] = save_np_datatype
    dispatch[numpy.float16] = save_np_datatype
    dispatch[numpy.float32] = save_np_datatype
    dispatch[numpy.float64] = save_np_datatype
    dispatch[numpy.complex_] = save_np_datatype
    dispatch[numpy.complex64] = save_np_datatype
    dispatch[numpy.complex128] = save_np_datatype


class DictToModel:
    """
    Rebuild a scikit-learn model from dict data generated by ModelToDict.save
    """

    def __init__(self):
        """ Store newly-built object
        """
        self.memo = {}

    def memoize(self, obj):
        l = len(self.memo)
        self.memo[l] = obj

    def load(self, data):
        return self.load_all(data[_OBJ])

    def load_all(self, data):
        """
        The main method to generate an object from dict data
        """
        dispatch = self.dispatch

        t = type(data)
        if t is dict:
            if _MEMO in data:
                return self.memo[data[_MEMO]]
            if _BYTES in data:
                return self.load_bytes(data[_BYTES])
            if _REDUCE in data:
                return self.load_reduce(data[_REDUCE])
            if _GLOBAL in data:
                return self.load_global(data[_GLOBAL])
            if _TUPLE in data:
                return self.load_tuple(data[_TUPLE])
            if _SET in data:
                return self.load_set(data[_SET])
            if _NP_NDARRAY in data:
                return self.load_np_ndarray(data[_NP_NDARRAY])
            if _NP_DATATYPE in data:
                return self.load_np_datatype(data[_NP_DATATYPE])
            return self.load_dict(data)
        f = dispatch.get(t)
        if f:
            return f(self, data)
        else:
            raise JsonPicklerError("Unsupported data found: %s" % str(data))

    dispatch = {}

    def load_primitive(self, data):
        return data

    dispatch[type(None)] = load_primitive
    dispatch[bool] = load_primitive
    dispatch[int] = load_primitive
    if six.PY2:
        dispatch[long] = load_primitive
    dispatch[float] = load_primitive
    dispatch[complex] = load_primitive

    def load_string(self, data):
        self.memoize(data)
        return data

    dispatch[str] = load_string

    def load_unicode(self, data):
        """
        `json.load` loads string as unicode in python 2,
        while some classes don't support unicode, like numpy.dtype
        """
        data = str(data)
        self.memoize(data)
        return data

    if six.PY2:
        dispatch[unicode] = load_unicode

    def load_bytes(self, data):
        data = data.encode('utf-8')
        self.memoize(data)
        return data

    def load_list(self, data):
        return [self.load_all(e) for e in data]

    dispatch[list] = load_list

    def load_tuple(self, data):
        obj = self.load_all(data)
        return tuple(obj)

    def load_set(self, data):
        obj = self.load_all(data)
        return set(obj)

    def load_dict(self, data):
        newdict = {}
        _keys = data[_KEYS]
        for k in _keys:
            try:
                v = data[k]
            # JSON dumps non-string key to string
            except KeyError:
                v = data[str(k)]
            newdict[k] = self.load_all(v)
        return newdict

    def load_global(self, data):
        module = data[0]
        name = data[1]
        func = SafePickler.find_class(module, name)
        self.memoize(func)
        return func

    def load_reduce(self, data):
        """
        Build object
        """
        _func = data[_FUNC]
        func = self.load_all(_func)
        assert callable(func), "%r" % func

        _args = data[_ARGS][_TUPLE]
        args = tuple(self.load_all(_args))

        try:
            obj = args[0].__new__(args[0], * args)
        except:
            obj = func(*args)

        _state = data.get(_STATE)
        if _state:
            state = self.load_all(_state)
            setstate = getattr(obj, "__setstate__", None)
            if setstate:
                setstate(state)
            else:
                assert (type(state) is dict)
                for k, v in state.items():
                    setattr(obj, k, v)

        self.memoize(obj)
        return obj

    def load_np_ndarray(self, data):
        _dtype = self.load_all(data[_DTYPE])
        _values = self.load_all(data[_VALUES])
        obj = numpy.array(_values, dtype=_dtype)
        return obj

    def load_np_datatype(self, data):
        _datatype = self.load_all(data[_DATATYPE])
        _value = self.load_all(data[_VALUE])
        obj = _datatype(_value)
        return obj


class H5Model:
    """
    Convert the dictionary model to h5 and reverse
    """

    def _create_dataset(self, file_obj, key, value):
        """
        Create dataset
        """
        try:
            file_obj.create_dataset(key, data=json.dumps(value))
        except:
            file_obj.create_dataset(key, data=value)

    def save_model(self, model, file_name):
        """
        Save the dictionary model to h5 file
        """
        serialised_model = ModelToDict().dump(model)
        h5file = h5py.File(file_name, 'w')

        # nested method for recursion
        def recursive_save_model(h5file_obj, dictionary):
            for model_key, model_value in dictionary.items():
                type_name = type(model_value)
                try:
                    if type_name is list:
                        if len(model_value) > 0:
                            list_obj = all(isinstance(x, dict) for x in model_value)
                            if list_obj is False:
                                self._create_dataset(h5file_obj, model_key, model_value)
                            else:
                                for index, model_item in enumerate(model_value):
                                    model_key_item = model_key + "/" + str(index)
                                    if model_item is not None:
                                        if model_key_item in h5file_obj:
                                            recursive_save_model(model_key_item, model_item)
                                        else:
                                            group = h5file_obj.create_group(model_key_item)
                                            recursive_save_model(group, model_item)
                                    else:
                                        self._create_dataset(h5file_obj, model_key_item, model_item)
                        else:
                            self._create_dataset(h5file_obj, model_key, model_value)
                    elif type_name is dict:
                        if isinstance(model_key, (int, numpy.int16, numpy.int32, numpy.int64)):
                            model_key = str(model_key)
                        if model_key in h5file_obj:
                            recursive_save_model(h5file_obj[model_key], model_value)
                        else:
                            group = h5file_obj.create_group(model_key)
                            recursive_save_model(group, model_value)
                    else:
                        self._create_dataset(h5file_obj, model_key, model_value)
                except Exception:
                    continue
        recursive_save_model(h5file, serialised_model)

    def _restore_list_in_model(self, model_object):
        """
        Convert the model dictionary to list if the keys are represented as numbers
        """
        if isinstance(model_object, dict):
            for key, value in model_object.items():
                if type(value) is dict:
                    keys = value.keys()
                    all_keys_number = all(str.isnumeric(x) for x in keys)
                    if all_keys_number is True:
                        keys = [int(ky) for ky in keys]
                        model_object[key] = list()
                        for idx in range(len(keys)):
                            model_object[key].append(value[str(idx)])
                self._restore_list_in_model(value)
        elif isinstance(model_object, list):
            for v in model_object:
                self._restore_list_in_model(v)

    def load_model(self, file_name):
        """
        Read the h5 file and reconstruct the model
        """
        reconstructed_model = dict()
        h5file = h5py.File(file_name, 'r')

        # nested method for recursion
        def recursive_load_model(h5file_obj, reconstructed_model):
            for key in h5file_obj.keys():
                # recurse if the item is a group
                if h5file_obj.get(key).__class__.__name__ is 'Group':
                    reconstructed_model[key] = dict()
                    recursive_load_model(h5file_obj[key], reconstructed_model[key])
                else:
                    try:
                        key_value = h5file_obj.get(key).value
                        reconstructed_model[key] = json.loads(key_value)
                    except Exception:
                        reconstructed_model[key] = key_value
                        continue
        recursive_load_model(h5file, reconstructed_model)
        H5Model()._restore_list_in_model(reconstructed_model)
        return DictToModel().load(reconstructed_model)
