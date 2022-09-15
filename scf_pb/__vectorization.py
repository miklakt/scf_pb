from typing import Callable
import inspect
import itertools
import numpy as np

def vectorize(product = True) -> Callable:
    def decorator(func : Callable) -> Callable:
        signature = inspect.signature(func)
        def wrapped(*args, **kwargs):
            args_to_kwargs = {param : arg for param, arg in zip(signature.parameters, args)}
            all_args = dict(**args_to_kwargs, **kwargs)
            vector_args = {k : v for k, v in all_args.items() if isinstance(v, (list, np.ndarray))}
            shape = [len(v) for v in vector_args.values()]
            if not shape:
                return func(**all_args)
            scalar_args = {k : v for k, v in all_args.items() if k not in vector_args}

            if product:
                iteration = itertools.product(*vector_args.values())
            #else:
            #    iteration = zip(*vector_args.values())

            results = np.empty(shape = shape, dtype = object)
            flatiter = results.flat
            for it in iteration:
                iteration_args = {k : v for k, v  in zip(vector_args, it)}
                iteration_args.update(scalar_args)
                results[flatiter.coords] = func(**iteration_args)
                next(flatiter)
                #results.append(result)
            #results = np.reshape(results, shape)
            #if not shape:
            #    results = results.item()
            return results
        wrapped.__signature__ = signature
        wrapped.__name__ = func.__name__
        wrapped.__module__ = func.__module__
        return wrapped
    return decorator

