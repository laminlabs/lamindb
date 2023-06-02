# %% [markdown]
# # Testing LaminDB's select statements
from datetime import datetime
from pathlib import Path
from typing import List, Optional

import sqlmodel as sqm

import lamindb as ln


def test_sqlmodel_select():
    if ln._USE_DJANGO:
        return None

    ln.setup.init(storage="sqlapi")

    class Hero(sqm.SQLModel, table=True):
        id: int = sqm.Field(primary_key=True)
        name: str
        team_id: Optional[int] = sqm.Field(foreign_key="team.id")
        team: "Team" = sqm.Relationship(back_populates="heroes")
        created_at: datetime = sqm.Field(default_factory=datetime.utcnow)

    class Team(sqm.SQLModel, table=True):
        id: int = sqm.Field(primary_key=True)
        name: str
        heroes: List[Hero] = sqm.Relationship(back_populates="team")

    Path("sqlapi/sqlapi.lndb").unlink()

    engine = sqm.create_engine("sqlite:///sqlapi/sqlapi.lndb")

    sqm.SQLModel.metadata.create_all(engine)

    with sqm.Session(engine) as session:
        team_1 = Team(id=0, name="stars")
        session.add(team_1)
        team_2 = Team(id=1, name="moons")
        session.add(team_2)
        session.add(Hero(id=0, name="test1", team_id=team_1.id))
        session.add(Hero(id=1, name="test2", team_id=team_1.id))
        session.add(Hero(id=2, name="test3"))
        session.commit()

    # ## Test select statements

    # %%
    ln.select(Hero, id=0).one()

    # %%
    ln.select(Hero).df()

    # %%
    ln.select(Hero, team_id=0).all()

    # %%
    ln.select(Hero).where(Hero.name == "test1").df()

    # %%
    ln.select(Hero).where(Hero.name == "test1").df()

    # %%
    ln.select(Hero).where(Hero.id == 0, Hero.name == "test1").df()

    # %%
    ln.select(Hero).where(Hero.id == 0, Hero.created_at <= datetime.utcnow()).df()

    # %%
    ln.select(Hero).where(sqm.or_(Hero.name == "test1", Hero.name == "test2")).df()

    # %%
    ln.select(Hero).where(sqm.or_(Hero.name == "test1", Hero.name == "test2")).offset(
        1
    ).df()

    # %%
    ln.select(Hero).where(sqm.or_(Hero.name == "test1", Hero.name == "test2")).offset(
        0
    ).limit(1).df()

    # %%
    ln.select(Hero).where(sqm.or_(Hero.name == "test1", Hero.name == "test2")).order_by(
        sqm.desc(Hero.created_at)
    ).df()

    # %%
    ln.select(Hero).join(Team).df()

    # %%
    ln.select(Hero, Team).join(Team).all()

    # %%
    ln.select(Hero, Team).join(Team).df()

    # %%
    ln.select(Hero, Team).where(Hero.team_id == Team.id).all()

    # %% [markdown]
    # ## Test autoflush

    # %%
    hero = Hero(name="test4")

    # %%
    hero

    # %%
    assert ln.select(Hero, name="test4").one_or_none() is None

    # %%
    team = ln.select(Team, name="stars").one()

    # %%
    hero.team = team

    # %%
    assert ln.select(Hero, name="test4").one_or_none() is None

    ln.setup.delete("sqlapi")
    ln.setup.load("lamindb-unit-tests")
