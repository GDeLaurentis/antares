import diskcache
import functools
import inspect
import hashlib
import pickle


class DiskCached(object):
    """DiskCache interface to Object."""

    CACHE_PATH = None
    DISKCACHE_SIZE_LIMIT_IN_GB = 32

    @property
    def diskcache(self):
        if self.CACHE_PATH is None:
            raise ValueError("CACHE_PATH is not set for DiskCached.")
        return diskcache.Cache(
            directory=self.CACHE_PATH,
            size_limit=self.DISKCACHE_SIZE_LIMIT_IN_GB * 2 ** 30
        )

    def summarize_diskcache(self, print_warnings=False):
        """Summarizes the status of the diskcache."""
        print(f"Cache at {self.diskcache.volume() / self.diskcache.size_limit * 100:.3}% of maximum size.")
        warnings = self.diskcache.check(fix=False)
        print(f"Number of warning: {len(warnings)}")  # if this is not zero there may be some corruption
        if print_warnings:
            for i, warning in enumerate(warnings):
                print(f"Warning #{i}:", warning)
        print(f"Number of phase-space points currently stored: {len(list(self.diskcache.iterkeys()))}")

    def memoized(*decorator_args, **decorator_kwargs):
        """Diskcaching decorator generator."""
        def memoized_decorator(func):
            signature = inspect.signature(func)
            default_kwargs = {
                kw: value.default for kw, value in signature.parameters.items() if value.default != inspect.Signature.empty
            }

            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                cached = kwargs.get("cached", False)
                verbose = kwargs.get("verbose", False)
                kwargs = default_kwargs | kwargs

                # Clean up kwargs and args based on 'ignore'
                cleaned_kwargs = {k: v for k, v in kwargs.items() if k not in decorator_kwargs.get('ignore', [])}
                ignore_args = decorator_kwargs.get('ignore', set())
                ignore_positions = {i for i in ignore_args if isinstance(i, int)}
                cleaned_args = [arg for i, arg in enumerate(args) if i not in ignore_positions]

                # Accessing default kwargs in here is non trivial, enche the inspect.signature story above
                if verbose:
                    print("decorator_args:", decorator_args)
                    print("decorator_kwargs:", decorator_kwargs)
                    print("args:", args)
                    print("default kwargs:", default_kwargs)
                    print("complete kwargs:", kwargs)
                    print("cleaned args:", cleaned_args)
                    print("cleaned kwargs:", cleaned_kwargs)

                if not cached:
                    if verbose:
                        print("Evaluating without cache")
                    return func(*args, **kwargs)

                self = args[0]

                if verbose:
                    print("Evaluating with cache")

                # Compute SHA-256 hash of the serialized cleaned args and kwargs
                cache_key_raw = (tuple(cleaned_args), frozenset(cleaned_kwargs.items()))
                cache_key_bytes = pickle.dumps(cache_key_raw)
                cache_key_hash = hashlib.sha256(cache_key_bytes).hexdigest()

                result = self.diskcache.get(cache_key_hash)
                if result is not None:
                    if verbose:
                        print("Read from cache")
                    return result

                if verbose:
                    print("Computing from scratch")
                result = func(*args, **kwargs)
                self.diskcache.add(cache_key_hash, result)
                return result

            return wrapper
        return memoized_decorator
