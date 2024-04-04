from __future__ import annotations


def fake_bio_notebook_titles(n=100) -> list[str]:
    """A fake collection of study titles."""
    from faker import Faker

    fake = Faker()

    from faker_biology.mol_biol import Antibody
    from faker_biology.physiology import CellType, Organ, Organelle

    fake.add_provider(CellType)
    fake.add_provider(Organ)
    fake.add_provider(Organelle)
    fake.add_provider(Antibody)

    my_words = [
        "study",
        "investigate",
        "research",
        "result",
        "cluster",
        "rank",
        "candidate",
        "visualize",
        "efficiency",
        "classify",
    ]
    my_words += [fake.organ() for i in range(5)] + ["intestine", "intestinal"]
    my_words += [fake.celltype() for i in range(10)]
    my_words += [fake.antibody_isotype() for i in range(20)]

    my_notebook_titles = [fake.sentence(ext_word_list=my_words) for i in range(n)]

    return my_notebook_titles
