"""Fields.

.. autosummary::
   :toctree: .

   CharField
   TextField
   ForeignKey
   BooleanField
   DateField
   DateTimeField
   BigIntegerField
   IntegerField
   OneToOneField
   FloatField
   DecimalField
   BinaryField
   JSONField
   EmailField
   TimeField
   SlugField
   URLField
   UUIDField
   PositiveIntegerField
   PositiveSmallIntegerField
   SmallIntegerField
   GenericIPAddressField
   DurationField
   CharField
   TextField
"""

from django.db import models


class CharField(models.CharField):
    """Custom `CharField` with default values for `blank`, `default`, and `max_length`.

    Django default values for `CharField` are `blank=False`, `default=""`, undefined `max_length`.
    """

    def __init__(self, max_length: int = 255, **kwargs):
        kwargs["max_length"] = max_length  # Set max_length in kwargs
        kwargs.setdefault("blank", True)
        kwargs.setdefault("default", None)
        super().__init__(**kwargs)  # Pass all arguments as kwargs


class TextField(models.TextField):
    """Custom `TextField` with default values for `blank` and `default`.

    Django default values for `TextField` are `blank=False`, `default=''`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        kwargs.setdefault("default", None)
        super().__init__(*args, **kwargs)


class ForeignKey(models.ForeignKey):
    """Custom `ForeignKey` with default values for `blank`.

    Django default value for `ForeignKey` `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


# fix doc string that otherwise errors
ForeignKey.get_extra_descriptor_filter.__doc__ = (
    ForeignKey.get_extra_descriptor_filter.__doc__.replace(
        ".filter(**kwargs)", "`.filter(**kwargs)`"
    )
)


class BooleanField(models.BooleanField):
    """Custom `BooleanField` with default values for `blank` and `default`.

    Django default values for `BooleanField` are `blank=False`, `default=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        kwargs.setdefault("default", None)
        super().__init__(*args, **kwargs)


class DateField(models.DateField):
    """Custom `DateField` with default values for `blank`.

    Django default values for `DateField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class DateTimeField(models.DateTimeField):
    """Custom `DateTimeField` with default values for `blank`.

    Django default values for `DateTimeField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class BigIntegerField(models.BigIntegerField):
    """Custom `BigIntegerField` with default values for `blank`.

    Django default values for `BigIntegerField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        kwargs.setdefault("default", None)
        super().__init__(*args, **kwargs)


class IntegerField(models.IntegerField):
    """Custom `IntegerField` with default values for `blank`.

    Django default values for `IntegerField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class OneToOneField(models.OneToOneField):
    """Custom `OneToOneField` with default values for `blank`.

    Django default values for `OneToOneField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class FloatField(models.FloatField):
    """Custom `FloatField` with default values for `blank`.

    Django default values for `FloatField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class DecimalField(models.DecimalField):
    """Custom `DecimalField` with default values for `blank`.

    Django default values for `DecimalField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class JSONField(models.JSONField):
    """Custom `JSONField` with default values for `blank`.

    Django default values for `JSONField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class DurationField(models.DurationField):
    """Custom `DurationField` with default values for `blank`.

    Django default values for `DurationField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class URLField(models.URLField):
    """Custom `URLField` with default values for `blank`.

    Django default values for `URLField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class EmailField(models.EmailField):
    """Custom `EmailField` with default values for `blank`.

    Django default values for `EmailField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class TimeField(models.TimeField):
    """Custom `TimeField` with default values for `blank`.

    Django default values for `TimeField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class SlugField(models.SlugField):
    """Custom `SlugField` with default values for `blank`.

    Django default values for `SlugField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class UUIDField(models.UUIDField):
    """Custom `UUIDField` with default values for `blank`.

    Django default values for `UUIDField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class PositiveIntegerField(models.PositiveIntegerField):
    """Custom `PositiveIntegerField` with default values for `blank`.

    Django default values for `PositiveIntegerField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class PositiveSmallIntegerField(models.PositiveSmallIntegerField):
    """Custom `PositiveSmallIntegerField` with default values for `blank`.

    Django default values for `PositiveSmallIntegerField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class SmallIntegerField(models.SmallIntegerField):
    """Custom `SmallIntegerField` with default values for `blank`.

    Django default values for `SmallIntegerField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class BinaryField(models.BinaryField):
    """Custom `BinaryField` with default values for `blank`.

    Django default values for `BinaryField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)


class GenericIPAddressField(models.GenericIPAddressField):
    """Custom `GenericIPAddressField` with default values for `blank`.

    Django default values for `GenericIPAddressField` are `blank=False`.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("blank", True)
        super().__init__(*args, **kwargs)
