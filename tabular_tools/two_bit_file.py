# class supporting reads from UCSC's .2bit format for DNA sequence collections.
#
# The .2bit file format spec can be found at
#   http://genome.ucsc.edu/FAQ/FAQformat#format7
#
# History: this was originally written circa 2005 by Bob Harris, and some
# version of it was incorporated in James Taylor's bx python package.

import sys

from struct      import *
from collections import MutableMapping
from sys         import stderr


class TwoBitFile(MutableMapping):

	TWOBIT_MAGIC_NUMBER      = 0x1A412743
	TWOBIT_MAGIC_NUMBER_SWAP = 0x4327411A
	TWOBIT_VERSION           = 0

	byteToDna = [None] * 256
	for ix in range(256):
		nuc1 = (ix >> 6)
		nuc2 = (ix >> 4) & 3
		nuc3 = (ix >> 2) & 3
		nuc4 =  ix       & 3
		byteToDna[ix] = ["TCAG"[nuc1],"TCAG"[nuc2],"TCAG"[nuc3],"TCAG"[nuc4]]

	def __init__(self,f,unmask=False,debug=None):
		self.debug = debug if (debug != None) else []

		# read magic and determine byte order

		if (type(f) == str):
			f = open(f,"rb")

		magic = unpack(">L",f.read(4))[0]
		if (magic == TwoBitFile.TWOBIT_MAGIC_NUMBER):
			self.byteOrder = ">"
		elif (magic == TwoBitFile.TWOBIT_MAGIC_NUMBER_SWAP):
			self.byteOrder = "<"
		else:
			raise Exception("Not a 2bit file")

		self.file   = f
		self.unmask = unmask

		self.version = self.read("L")
		if (self.version != TwoBitFile.TWOBIT_VERSION):
			raise Exception("unsupported 2bit version %d" % self.version)

		# read the sequence names and create a hash

		self.seqCount = self.read("L")
		reserved      = self.read("L")
		if (reserved != 0):
			raise Exception("2bit file has %d at offset 0x0000000C (supposed to be zero)" \
			              % reserved)
		readOffset = 0x0000000C + 4

		self.index = {}
		for nameIx in range(self.seqCount):
			nameOffset = readOffset
			name = self.read_p_string()
			headerOffset = self.read("L")
			readOffset += 1 + len(name) + 4
			if (name in self.index):
				raise Exception("sequence name \"%s\" appears more than once in 2bit file index (fileIndex[%d] at file offset 0x%08X)" \
				              % (name,nameIx,nameOffset))
			self.index[name] = headerOffset  # this will be replaced by seq object if/when fetched
		self.preloadedCount = 0

		self.len = {}
		for name in self.index:
			self.len[name] = self.sequence_length(name)

	def __getitem__(self,name):
		try:
			seq = self.index[name]
			if (type(seq) == int):
				self.index[name] = seq = TwoBitSequence(self,seq)
			if (seq.ready): return seq
			self.preload_sequence(name)
			return seq
		except KeyError:
			return None

	def __setitem__(self,name,val):
		raise Exception("TwoBitFile doesn't support __setitem__)")

	def __delitem__(self,name):
		raise Exception("TwoBitFile doesn't support __delitem__)")

	def __len__(self):
		return len(self.index)

	def __iter__(self):
		for name in self.index:
			yield name

	def keys(self):
		return self.index.keys()

	def sequence_length(self,name):
		seq = self.index[name]
		if (type(seq) != int):
			# the sequence has been preloaded; ask it for its length
			return len(seq)
		else:
			# the sequence has not been preloaded; read its length from file
			headerOffset = seq
			self.file.seek(headerOffset)
			return self.read("L")

	def preload_sequence(self,name):
		seq = self.index[name]
		if (type(seq) == int):
			self.index[name] = seq = TwoBitSequence(self,seq)

		self.file.seek(seq.headerOffset)
		seq.size = self.read("L")
		if ("blocks" in self.debug):
			print ("reading blocks for \"%s\" (%d/%d)" \
			     % (name,self.preloadedCount+1,self.seqCount),
			       file=stderr)
		nBlocks    = self.read_blocks_info()
		maskBlocks = self.read_blocks_info()
		reserved   = self.read("L")
		if (reserved != 0):
			raise Exception("2bit file has %d at offset 0x%08X (supposed to be zero)" \
			              % (reserved,self.file.tell()-4))
		seq.sequenceOffset = self.file.tell()

		seq.nBlocks = nBlocks
		if (self.unmask): seq.maskBlocks = []
		else:             seq.maskBlocks = maskBlocks

		seq.ready = True
		self.preloadedCount += 1

	def read_blocks_info(self):
		blockCount = self.read("L")
		if (blockCount == 0): return []
		starts = self.read(str(blockCount)+"L",untuple=False)
		sizes  = self.read(str(blockCount)+"L",untuple=False)
		blocks = [(start,start+size) for (start,size) in zip(starts,sizes)]
		blocks.sort()
		return blocks

	def read(self,pattern,untuple=True):
		pattern = self.byteOrder + pattern
		val = unpack(pattern,self.file.read(calcsize(pattern)))
		if (untuple) and (len(val) == 1): return val[0]
		else:                             return val

	def read_p_string(self):
		length = self.read("B")
		return self.file.read(length).decode("utf-8")


class TwoBitSequence(object):
	def __init__(self,tbf,headerOffset):
		self.tbf            = tbf
		self.headerOffset   = headerOffset
		self.size           = None
		self.nBlocks        = None  # sorted list of (start,end)
		self.maskBlocks     = None  # sorted list of (start,end)
		self.sequenceOffset = None
		self.ready          = False

	def __getitem__(self,slice):
		if ("getitem" in self.tbf.debug):
			print ("slice = %s" % slice,file=stderr)
		(start,stop,stride) = slice.indices(self.size)
		assert (stride == 1), "Striding in slices not supported"
		if (stop <= start): return ""
		return self.read_slice(start,stop)

	def __len__(self):
		return self.size

	def get(self,start,end):
		if (start < 0):         start = 0
		if (end   > self.size): end  = self.size
		if (stop <= start): return ""
		return self.read_slice(start,stop)

	def read_slice(self,start,end):
		# $$$ we could improve this by checking the nBlocks first to see if we
		# $$$ .. will be N-masking the whole slice;  save the relevant part of
		# $$$ .. that list (if any), read, then apply

		if ("slice" in self.tbf.debug):
			print("start,end = %s,%s" % (start,end),file=stderr)

		nucsToFetch = end - start

		packedStart = start >> 2
		packedEnd   = (end+3) >> 2
		packedSize  = packedEnd - packedStart

		self.tbf.file.seek(self.sequenceOffset + packedStart)
		packedBytes = self.tbf.file.read(packedSize)

		if (packedSize == 1):
			# dna is entirely within one packed byte
			sOff = start - (packedStart << 2)
			eOff = sOff  + nucsToFetch
			packed   = packedBytes[0]
			unpacked = TwoBitFile.byteToDna[packed]
			dna = unpacked[sOff:eOff]
		else:
			# dna spans many packed bytes
			dna = []
			packedIx = 0

			skip = start & 3
			if (skip > 0):
				keep = 4 - skip
				packed   = packedBytes[packedIx]
				unpacked = TwoBitFile.byteToDna[packed]
				packedIx += 1
				dna += unpacked[skip:]
				nucsToFetch -= keep

			while (nucsToFetch >= 4):
				packed = packedBytes[packedIx]
				dna    += TwoBitFile.byteToDna[packed]
				packedIx += 1
				nucsToFetch -= 4

			if (nucsToFetch > 0):
				packed   = packedBytes[packedIx]
				unpacked = TwoBitFile.byteToDna[packed]
				dna += unpacked[:nucsToFetch]

		# apply N- and mask-intervals

		for (s,e) in self.nBlocks:
			if (e <= start): continue
			if (s >= end):   break
			if (s < start):  s = start
			if (e > end):    e = end
			s -= start
			e -= start
			dna[s:e] = ["N"] * (e-s)

		for (s,e) in self.maskBlocks:
			if (e <= start): continue
			if (s >= end):   break
			if (s < start):  s = start
			if (e > end):    e = end
			s -= start
			e -= start
			dna[s:e] = [ch for ch in "".join(dna[s:e]).lower()]

		return "".join(dna)

