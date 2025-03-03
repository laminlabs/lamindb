from typing import TYPE_CHECKING, Literal, Union, get_args, get_origin, get_type_hints

from lamin_utils import colors

from lamindb.errors import FieldValidationError

if TYPE_CHECKING:
    from ..models import Record


def validate_literal_fields(record: "Record", kwargs) -> None:
    """Validate all Literal type fields in a record.

    Args:
        record: record being validated

    Raises:
        ValidationError: If any field value is not in its Literal's allowed values
    """
    # check is based on string to avoid circular imports
    if record.__class__.__name__ == "Feature":
        # the FeatureDtype is more complicated than a simple literal
        # because it allows constructs like cat[ULabel] etc.
        # the User model is used at startup and throws a datetime-related error otherwise
        # simmilar for Storage & Source
        return None
    try:
        type_hints = get_type_hints(record.__class__)
    except TypeError:
        # for 3.9, get_type_hints errors with | in type hints
        return
    errors = {}

    for field_name, field_type in type_hints.items():
        # Handle both plain Literal and Union/Optional Literal types
        origin = get_origin(field_type)
        if origin is Union:
            # For Optional/Union types, find the Literal type if it exists
            literal_type = next(
                (t for t in get_args(field_type) if get_origin(t) is Literal), None
            )
        else:
            # For plain types, check if it's a Literal
            literal_type = field_type if origin is Literal else None

        # Skip if no Literal type found
        if literal_type is None:
            continue

        value = kwargs.get(field_name)
        if value is not None:
            valid_values = set(get_args(literal_type))
            if value not in valid_values:
                errors[field_name] = (
                    f"{field_name}: {colors.yellow(value)} is not a valid value"
                    f"\n    → Valid values are: {colors.green(', '.join(sorted(valid_values)))}"
                )

    if errors:
        message = "\n  "
        for _, error in errors.items():
            message += error + "\n  "
        raise FieldValidationError(message)
