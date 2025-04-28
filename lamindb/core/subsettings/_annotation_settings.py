class AnnotationSettings:
    n_max_records: int = 1000
    """Maximal number of records to annotate with during automated annotation.

    If the number of records to annotate exceeds this limit, print a warning and do not annotate.

    The number is calculated per feature for labels, and per schema for features.
    """


annotation_settings = AnnotationSettings()
