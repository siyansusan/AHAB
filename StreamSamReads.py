# =========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
# Author: Siyan Liu
# =========================================================================
from __future__ import (absolute_import, division, print_function,
                        unicode_literals, generators, nested_scopes, with_statement)

import logging

from SamReader import SamReader
from SamReadGroup import SamReadGroup


class StreamSamReads:
    """
    This is an adapter class that reads SamRecords from a SamReader and groups
    them into SamPairedRead objects.  It implements a buffer, to avoid losing
    reads when reading too far into the SAM file.

    Attributes:
        reader : SamReader
        dedup : boolean
        bufferedRec : SamRecord
        bufferedPair : SamPairedRead
    Instance Methods:
        stream=SamPairedReadStream(filename,dedup=True)
        pair=stream.nextPair() # returns SamPairedRead
        readGroup=stream.nextGroup() # returns array of SamPairedRead
    Class Methods:
        none
    """

    def __init__(self, filename, dedup=True):
        self.reader = SamReader(filename)
        self.dedup = dedup
        self.buffer_read = None

    def nextGroup(self):
        group = SamReadGroup()
        readID = None
        while True:

            # Obtain next read sequence
            if self.buffer_read == None:
                read = self.reader.nextSequence()
            else:
                read = self.buffer_read
                self.buffer_read = None

            # Check if any read was read
            if read is None:
                break

            # Skip if read is unmapped
            if read.flag_unmapped():
                logging.debug("Read is unmapped")
                continue

            # Skip if read is marked as a PCR duplicate
            if self.dedup and read.flag_PCRduplicate():
                logging.debug("Read is PCR duplicate")
                continue

            # Record the first read ID from the group
            if readID is None:
                readID = read.getID()

            # Add read to group if it has the same ID as the first one
            if read.getID() == readID:
                group.reads.append(read)
            else:
                self.buffer_read = read
                break

        return group
