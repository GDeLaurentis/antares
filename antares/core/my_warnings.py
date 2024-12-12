from __future__ import print_function
from __future__ import unicode_literals

import logging
import multiprocessing


class warnings(object):

    def __init__(self):
        self._warning = None
        self.manager = multiprocessing.Manager()
        self._silent = self.manager.Value('b', False)
        self.past = self.manager.list()

    def warn(self, message, warning=None):
        self._warning = warning
        # log the warning only if it is not the same as the previous one (avoid spam)
        if self.previous_message != message:
            logging.warning(message)
        self.previous_message = message
        if self.silent is False and self._warning not in self.past:
            print("\nWarning:" + message)
        if warning is not None and warning not in self.past:
            self.past += [warning]

    def raise_error(self, message, warning=None):
        self._warning = warning
        logging.error(message)
        print("\nError:" + message)
        if warning is not None and warning not in self.past:
            self.past += [warning]

    def clear(self):
        self._warning = None

    @property
    def warning(self):
        return self._warning

    @property
    def silent(self):
        return self._silent.value

    @silent.setter
    def silent(self, temp_value):
        self._silent.value = temp_value


oWarning = warnings()
oWarning.previous_message = ""
