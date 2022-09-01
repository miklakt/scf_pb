from scf_pb import _scf_pb
from scf_pb.__vectorization import vectorize
import inspect
from typing import Any, Callable, List

def lexical_cast_builtins(typename):
    try:
        type_ = getattr(__builtins__, typename)
    except AttributeError:
        type_ = typename
    return type_

def parse_signature(text_signature):

    func_name, rest = text_signature.split("(",1)
    args, return_type = rest.split(")", 1)
    return_type = return_type.rsplit(" ", 1)[-1]
    return_type = lexical_cast_builtins(return_type)

    args = [[item.strip() for item in arg.split(",")] for arg in args.split("*")]
    [item.remove('') for item in args]
    args = [[[it.strip() for it in item.split(':')] for item in arg] for arg in args]

    positional_or_keyword, keyword_only = [{k: lexical_cast_builtins(v) for k,v in arg} for arg in args]

    return func_name, return_type, positional_or_keyword, keyword_only

def signature_from_pybind_doc(doc):
    from inspect import Parameter, Signature
    text_sig = doc.split("\n")[0]
    func_name, return_type, positional_or_keyword, keyword_only = parse_signature(text_sig)
    positional_or_keyword = [Parameter(k, Parameter.POSITIONAL_OR_KEYWORD, annotation=v) for k,v in positional_or_keyword.items()]
    keyword_only = [Parameter(k, Parameter.KEYWORD_ONLY, annotation=v) for k,v in keyword_only.items()]

    sig = Signature([*positional_or_keyword, *keyword_only], return_annotation=return_type)
    return sig, func_name

def wrap_pybind_func(pybind_func):
    import functools
    doc = pybind_func.__doc__
    sig, func_name = signature_from_pybind_doc(doc)
    @functools.wraps(pybind_func)
    def wrapper(*args, **kwargs):
        return pybind_func(*args, **kwargs)
    wrapper.__signature__ = sig
    return wrapper

D_eff_cxx = wrap_pybind_func(_scf_pb.D_eff)
D_eff = vectorize()(D_eff_cxx)

    
    