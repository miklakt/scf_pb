from typing import Callable
import inspect
import itertools
import numpy as np
import os

import scf_pb

from concurrent.futures import ThreadPoolExecutor, as_completed

def eval_expressions_in_arguments(kwargs : dict):
    #make it lazy
    computables = {
        "H" : (lambda : scf_pb.D(**{ k:kwargs[k] for k in ["N", "sigma", "chi"]})),
        "D" : (lambda : scf_pb.D(**{ k:kwargs[k] for k in ["N", "sigma", "chi"]})),
        "phi" : (lambda : scf_pb.phi(**{ k:kwargs[k] for k in ["N", "sigma", "chi", "z"]})),
        }

    #add constants to locals
    constants = {k:v for k, v in kwargs.items() if not isinstance(v, str)}

    #one variable substitutuions
    singles_substitutions = {k:constants[v] for k, v in kwargs.items() if v in constants.keys()}
    constants.update(singles_substitutions)

    #the rest is expressions
    expressions = {k:v for k, v in kwargs.items() if ((isinstance(v, str)) and (v not in constants.keys()))}
    if expressions:
        computed = {}
        for k, expr in expressions.items():
            code = compile(expr, "<string>", "eval")
            for name in code.co_names:
                if name in computables.keys(): computed[name] = computables[name]()
            constants[k] = eval(code, {**constants, **computed})


    return constants



def vectorize() -> Callable:
    def decorator(func : Callable) -> Callable:
        signature = inspect.signature(func)

        def wrapped(*args, **kwargs):

            progressbar = kwargs.pop("progressbar", False)
            max_workers = kwargs.pop("max_workers", os.cpu_count())

            args_to_kwargs = {param : arg for param, arg in zip(signature.parameters, args)}
            all_args = dict(**args_to_kwargs, **kwargs)

            vector_args = {k : v for k, v in all_args.items() if isinstance(v, (list, np.ndarray))}
            shape = [len(v) for v in vector_args.values()]
            if not shape:
                return func(**all_args)
            scalar_args = {k : v for k, v in all_args.items() if k not in vector_args}

            iteration = itertools.product(*vector_args.values())

            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                futures = []
                for it in iteration:
                    iteration_args = {k : v for k, v  in zip(vector_args, it)}
                    iteration_args.update(scalar_args)
                    iteration_args = eval_expressions_in_arguments(iteration_args)
                    futures.append(pool.submit(func, **iteration_args))

                if progressbar:
                    import tqdm
                    pbar = tqdm.tqdm(total = np.product(shape), desc = f"{func.__name__}, {max_workers} tread(s)")
                    for future_completed in as_completed(futures):
                        pbar.update()

            results = np.empty(shape = shape, dtype = object)
            flatiter = results.flat
            for future in futures:
                results[flatiter.coords] = future.result()
                next(flatiter)

            return results

        wrapped.__signature__ = signature
        wrapped.__name__ = func.__name__
        wrapped.__module__ = func.__module__
        return wrapped
    return decorator
