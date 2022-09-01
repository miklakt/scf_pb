from typing import Callable
import inspect
import itertools
try:
    import numpy as np
except:
    pass

def vectorize(product = True, signature = None) -> Callable:
    def decorator(func : Callable) -> Callable:
        nonlocal signature
        if signature is None:
            try:
                signature = inspect.signature(func)
            except:
                raise TypeError("Can not find signature for the callable, please provide one")
        def wrapped(*args, **kwargs):
            args_to_kwargs = {param : arg for param, arg in zip(signature.parameters, args)}
            all_args = dict(**args_to_kwargs, **kwargs)
            vector_args = {k : v for k, v in all_args.items() if isinstance(v, (list, np.ndarray))}
            shape = [len(v) for v in vector_args.values()]
            scalar_args = {k : v for k, v in all_args.items() if k not in vector_args}

            if product:
                iteration = itertools.product(*vector_args.values())
            else:
                iteration = zip(*vector_args.values())

            results = []
            for it in iteration:
                iteration_args = {k : v for k, v  in zip(vector_args, it)}
                iteration_args.update(scalar_args)
                results.append(func(**iteration_args))
            results = np.reshape(results, shape)
            if not shape:
                results = results.item()
            return results
        wrapped.__signature__ = signature
        wrapped.__name__ = func.__name__
        wrapped.__module__ = func.__module__
        return wrapped
    return decorator