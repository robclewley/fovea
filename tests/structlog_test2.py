from __future__ import absolute_import, division, print_function

import structlog


class ListLogger(object):
    def __init__(self):
        self.events = []

    def add_event(self, **kw):
        self.events.append(kw)


logger = ListLogger()

log = structlog.wrap_logger(
    logger,
    processors=[structlog.processors.TimeStamper(fmt="iso")]
)

log.add_event("ford", towel=True)
log.add_event("trillian", towel=False)

print(logger.events)
