from __future__ import annotations

import lamindb as ln


def test_flow_pep604_union_annotation_keeps_valid_param():
    run_str = None
    run_list = None
    transform = None
    try:

        @ln.flow(global_run="clear")
        def typed_flow(parents: str | list[str]) -> str:
            assert ln.context.run is not None
            return ln.context.run.uid

        run_str = ln.Run.get(uid=typed_flow("root"))
        run_list = ln.Run.get(uid=typed_flow(["root", "child"]))
        transform = run_str.transform
        assert run_str.params == {"parents": "root"}
        assert run_list.params == {"parents": ["root", "child"]}
    finally:
        ln.context._run = None
        if run_str is not None:
            run_str.delete(permanent=True)
        if run_list is not None:
            run_list.delete(permanent=True)
        if transform is not None:
            transform.delete(permanent=True)


def test_flow_pep604_union_annotation_skips_invalid_param():
    run = None
    transform = None
    try:

        @ln.flow(global_run="clear")
        def typed_flow(parents: str | list[str]) -> str:
            assert ln.context.run is not None
            return ln.context.run.uid

        run = ln.Run.get(uid=typed_flow(42))
        transform = run.transform
        assert run.params == {}
    finally:
        ln.context._run = None
        if run is not None:
            run.delete(permanent=True)
        if transform is not None:
            transform.delete(permanent=True)


def test_flow_optional_annotation_from_future_annotations():
    run_none = None
    run_str = None
    transform = None
    try:

        @ln.flow(global_run="clear")
        def typed_flow(parents: str | None = None) -> str:
            assert ln.context.run is not None
            return ln.context.run.uid

        run_none = ln.Run.get(uid=typed_flow())
        run_str = ln.Run.get(uid=typed_flow("root"))
        transform = run_str.transform
        # None-valued params are intentionally omitted from run.params.
        assert run_none.params == {}
        assert run_str.params == {"parents": "root"}
    finally:
        ln.context._run = None
        if run_none is not None:
            run_none.delete(permanent=True)
        if run_str is not None:
            run_str.delete(permanent=True)
        if transform is not None:
            transform.delete(permanent=True)
