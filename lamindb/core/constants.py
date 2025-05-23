from enum import Enum


class Colors(str, Enum):
    """Main colors used throughout our stack.

    Colors that are prefixed with `LAMIN` are our branding colors.
    """

    LAMIN_GREEN = "#12b981"
    LAMIN_ORANGE = "#fb923c"
    LAMIN_BLUE = "#2F5F9E"
    GREEN_LIGHTER = "#10b981"
    GREEN_DARKER = "#065f46"
    GREEN_FILL = "honeydew"

    def __str__(self):
        return self.value
