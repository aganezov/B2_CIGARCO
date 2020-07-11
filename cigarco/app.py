import argparse
import sys
from argparse import Namespace
from dataclasses import dataclass, field, InitVar
import logging
from logging import Logger
from typing import Iterator, Optional, Tuple, Iterable

from cigarco.mapping import Alignment, TransformationQuery, CManager, TransformedResult


@dataclass
class AlignmentsStreamer(object):
    """
    Iterator-like object that provides generative behaviour for creating Alignment object from string representation, while ensuring data validation

    Args:
        alignments_str_source (Iterable[src]): a source of alignments represented as strings

    """
    alignments_str_source: Iterable[str]

    def __iter__(self):
        for str_entry in self.alignments_str_source:
            data = str_entry.strip().split("\t")
            if len(data) < 4:
                raise ValueError(f"Insufficient info alignment definition '{str_entry.strip()}'")
            q_name = data[0]
            t_name = data[1]
            try:
                coordinate = int(data[2])
                cigar = data[3]
                yield Alignment(q_name, t_name, coordinate, cigar)
            except ValueError:
                raise ValueError(f"Could not parse alignment definition '{str_entry.strip()}'")


@dataclass
class QueryStreamer(object):
    """
    Iterator-like object that provides generative behaviour for creating TransformationQuery objects from string representation, while ensuring data validation

    Args:
        alignments_str_source (Iterable[src]): a source of alignments represented as strings
    """
    queries_str_source: Iterable[str]

    def __iter__(self) -> Iterable[TransformationQuery]:
        for str_entry in self.queries_str_source:
            data = str_entry.strip().split("\t")
            if len(data) < 2:
                raise ValueError(f"Insufficient info in query definition '{str_entry.strip()}'")
            q_name = data[0]
            try:
                coordinate = int(data[1])
                yield TransformationQuery(q_name, coordinate)
            except ValueError:
                raise ValueError(f"Could not parse alignment definition '{str_entry.strip()}'")


@dataclass
class CigarcoApp(object):
    """ Application like class an instance of which is designed to take transformation input (both alignments and queries)
        and provides a generative function for the actual coordinate transformation functionality

        Includes logging facilities, though debug logging is not yet implemented
        Provides three levels (I, R, and F) for handling errors (both in the input alignments and in the coordinate transformation queries),

    Args:
        alignments (Iterable[Alignment]): and iterable collection of Alignment class instances based on which the coordinate transformation queries are to be evaluated
        queries (Iterable[Alignment]): and iterable collection of TransformationQuery class instance that determine an specific coordinate trasnforamtion task
        fail_mode (choice of ['I', 'R', 'F']): a single-char enum that specifies hto the app shall behave in an event of an encountered error
            (ignore, report in log, and fail, respectively)
    """
    alignments: Iterable[Alignment]
    queries: Iterable[TransformationQuery]
    _cmanager: CManager = field(init=False, repr=False)
    fail_mode: str = 'R'
    _logger: Logger = None
    log_level: InitVar[int] = logging.INFO

    @staticmethod
    def _setup_logger(logger: logging.Logger, level: int):
        """ A helper staticmethod designed to encapsulate logger setup
        """
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(asctime)s %(name)s [%(levelname)-8.8s] -- %(message)s")
        handler.setFormatter(formatter)
        handler.setLevel(level)
        logger.addHandler(handler)

    def __post_init__(self, log_level):
        """ Logic of logger setup and actually adding alignments to internal CManager instance are implemented here, because of the @dataclass nature of the Application class

        If errors are encountered during alignment iteration/addition, a fail_mode app attribute dictates if and hwo the app shall proceed
        """
        self._logger = logging.getLogger('CIGARCO')
        self._setup_logger(self._logger, log_level)
        self._logger.info("Starting CigarcoApp")
        self._logger.info("Adding alignment records to the mapping manager")
        self._cmanager = CManager()
        alignments: Iterator[Alignment] = iter(self.alignments)
        while True:
            try:
                alignment = next(alignments, None)
                if alignment is None:
                    break
                self._cmanager.add_alignment(alignment)
            except ValueError as e:
                if self.fail_mode == "I":
                    continue
                if self.fail_mode in ["R", "F"]:
                    self._logger.error(str(e))
                    if self.fail_mode == "F":
                        self._logger.critical("Exiting because of the error in input alignment data processing and fail mode set to 'F'. "
                                              "To allow for errors in the future to be skipped/reported use --error-mode 'I' or 'R'")
                        exit(1)

    def transformations_iter(self) -> Iterable[Tuple[TransformationQuery, TransformedResult]]:
        """ Generator-like wrapper for the coordinate transformation queries being actually execute
        All computational logic and matching of queries to alignments is outsourced to the the internal CManager instance
        Course of action in an even of encountered errors is determined by the `fail_mode` flag

        Returns:
            A tuple of TransformationQuery and the respective TransformationResult objects
        """
        queries: Iterator[TransformationQuery] = iter(self.queries)
        self._logger.info("Processing coordinate transformation queries")
        while True:
            try:
                query: Optional[TransformationQuery] = next(queries, None)
                if query is None:
                    break
                yield query, self._cmanager.transform_coordinate(query.seq_name, query.coordinate)
            except ValueError as e:
                if self.fail_mode == "I":
                    continue
                if self.fail_mode in ["R", "F"]:
                    self._logger.error(str(e))
                    if self.fail_mode == "F":
                        self._logger.critical("Exiting because of the error in coordinate transformation data processing and fail mode set to 'F'. "
                                              "To allow for errors in the future to be skipped/reported use --error-mode 'I' or 'R'")
                        exit(1)


def create_cli_parser():
    """
    Basic command line parser setup for using CIGARCO APp in console mode
    Returns:
        parser (argparse.ArgumentParser)
    """
    result = argparse.ArgumentParser(prog="CIGARCO",
                                     description="A utility application to transform position coordinates from query to target of alignments via a CIGAR string relationship", )
    result.add_argument("-a", "--alignments", type=argparse.FileType("rt"), required=True,
                        help="A tab-separated file with alignment records (1 per line) with a scheme 'Qname Tname start CIGAR'")
    result.add_argument("-q", "--queries", type=argparse.FileType("rt"), required=True,
                        help="A tab-separated file with transformation queries (1 per line) with a scheme 'Qname coordinate'")
    result.add_argument("--error-mode", choices=['I', 'R', 'F'], default='R',
                        help="Error mode for the application, where with 'I' errors are ignore/skipped, "
                             "with 'R' error are reported in the output, and with 'F' app crashes on first error; Default = R")
    result.add_argument("--log-level", choices=[0, 10, 20, 30, 40, 50], type=int, default=20)
    result.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout,
                        help="File to write the results of the coordinate transformations to. Default = stdout ")
    return result


def main(args_list):
    """ main entry point for the application run in a command line mode
    Args:
        args_list: sys.argv[1:] elements
    """
    parser = create_cli_parser()
    args: Namespace = parser.parse_args(args_list)
    alignment_streamer: AlignmentsStreamer = AlignmentsStreamer(args.alignments)
    query_streamer: QueryStreamer = QueryStreamer(args.queries)
    app = CigarcoApp(alignment_streamer, query_streamer, args.error_mode, log_level=args.log_level)
    for query, result in app.transformations_iter():
        print(query.seq_name, query.coordinate, result.seq_name, result.coordinate, sep="\t", file=args.output)


def execute_script():
    """
    wrapper around the main entry point for testing purposes
    """
    if __name__ == "__main__":
        main(sys.argv[1:])


execute_script()
