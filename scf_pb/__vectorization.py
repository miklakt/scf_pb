from typing import Callable
import inspect
import itertools
import numpy as np

from concurrent.futures import ThreadPoolExecutor

def vectorize() -> Callable:
    def decorator(func : Callable) -> Callable:
        signature = inspect.signature(func)
        def wrapped(*args, **kwargs):

            progressbar = kwargs.pop("progressbar", False)
            max_workers = kwargs.pop("max_workers", 1)

            args_to_kwargs = {param : arg for param, arg in zip(signature.parameters, args)}
            all_args = dict(**args_to_kwargs, **kwargs)
            vector_args = {k : v for k, v in all_args.items() if isinstance(v, (list, np.ndarray))}
            shape = [len(v) for v in vector_args.values()]
            if not shape:
                return func(**all_args)
            scalar_args = {k : v for k, v in all_args.items() if k not in vector_args}

            iteration = itertools.product(*vector_args.values())

            pool_results = []
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                for it in iteration:
                    iteration_args = {k : v for k, v  in zip(vector_args, it)}
                    iteration_args.update(scalar_args)
                    #results[flatiter.coords] = func(**iteration_args)
                    pool_results.append(pool.submit(func, **iteration_args))
                    #next(flatiter)
            
            results = np.empty(shape = shape, dtype = object)
            flatiter = results.flat
            
            if progressbar:
                import tqdm
                pool_results = tqdm.tqdm(pool_results)
            
            for pool_result in pool_results:
                results[flatiter.coords] = pool_result.result()
                next(flatiter)
            
            return results

        wrapped.__signature__ = signature
        wrapped.__name__ = func.__name__
        wrapped.__module__ = func.__module__
        return wrapped
    return decorator
