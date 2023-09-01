def attach_func_to_class_method(func_name, cls, globals):
    implementation = globals[func_name]
    target = getattr(cls, func_name)
    # assigning the original class definition docstring
    # to the implementation only has an effect for regular methods
    # not for class methods
    # this is why we need @doc_args for class methods
    implementation.__doc__ = target.__doc__
    setattr(cls, func_name, implementation)
