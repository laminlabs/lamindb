import lamindb as ln


def track_if_not_tracked(lamin_instance: str):
    ln.connect(lamin_instance)
    if not ln.context.run:
        ln.track()
    print(ln.context.run)
