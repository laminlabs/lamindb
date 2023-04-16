# Parts of this class are from the Scanpy equivalent, see license below

# BSD 3-Clause License

# Copyright (c) 2017 F. Alexander Wolf, P. Angerer, Theis Lab
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from contextlib import contextmanager
from enum import IntEnum
from logging import getLevelName
from typing import Any, ContextManager, Tuple, Union

from lamin_logger import RootLogger
from lamin_logger._logger import INFO, VERBOSITY_TO_LOGLEVEL, set_log_level


class Verbosity(IntEnum):
    error = 0
    warn = 1
    info = 2
    hint = 3
    debug = 4

    @property
    def level(self) -> int:
        # getLevelName(str) returns the int level…
        return getLevelName(VERBOSITY_TO_LOGLEVEL[self])

    @contextmanager  # type: ignore
    def override(self, verbosity: "Verbosity") -> ContextManager["Verbosity"]:  # type: ignore  # noqa
        """Temporarily override verbosity."""
        settings.verbosity = verbosity
        yield self
        settings.verbosity = self


def type_check(var: Any, varname: str, types: Union[type, Tuple[type, ...]]):
    if isinstance(var, types):
        return
    if isinstance(types, type):
        possible_types_str = types.__name__
    else:
        type_names = [t.__name__ for t in types]
        possible_types_str = "{} or {}".format(
            ", ".join(type_names[:-1]), type_names[-1]
        )
    raise TypeError(f"{varname} must be of type {possible_types_str}")


class _Settings:
    """Settings.

    For setup-related settings, see :class:`lamindb.setup.settings`.
    """

    def __init__(self):
        self._root_logger = RootLogger(INFO, self)

    error_on_file_hash_exists: bool = True
    """Upon ingestion, error if a file hash equals an existing hash in the DB.

    FAQ: :doc:`/faq/ingest-same-file-twice`.
    """
    track_run_inputs_upon_load: bool = False
    """Upon load, add loaded files as the input of the current notebook run.

    FAQ: :doc:`/faq/track-runin`.
    """

    @property
    def verbosity(self) -> Verbosity:
        """Verbosity level (default `warning`).

        Level 0: only show 'error' messages.
        Level 1: also show 'warning' messages.
        Level 2: also show 'info' messages.
        Level 3: also show 'hint' messages.
        Level 4: also show very detailed progress for debugging.
        """
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity: Union[Verbosity, int, str]):
        verbosity_str_options = [v for v in VERBOSITY_TO_LOGLEVEL if isinstance(v, str)]
        if isinstance(verbosity, Verbosity):
            self._verbosity = verbosity
        elif isinstance(verbosity, int):
            self._verbosity = Verbosity(verbosity)
        elif isinstance(verbosity, str):
            verbosity = verbosity.lower()
            if verbosity not in verbosity_str_options:
                raise ValueError(
                    f"Cannot set verbosity to {verbosity}. "
                    f"Accepted string values are: {verbosity_str_options}"
                )
            else:
                self._verbosity = Verbosity(verbosity_str_options.index(verbosity))
        else:
            type_check(verbosity, "verbosity", (str, int))
        set_log_level(self, VERBOSITY_TO_LOGLEVEL[self._verbosity])


settings = _Settings()
