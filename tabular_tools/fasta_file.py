"""
Note: the initial implementation of this predates pyfaidx. I maintain this to
      support backward compatibility for old scripts. But new scripts will do
      much better by using pyfaidx instead of this.
"""

from collections import MutableMapping
from os.path     import exists
from sys         import stderr

try:
	from hashlib import md5 as md5_new
except ImportError:
	try:
		from md5 import new as md5_new
	except ImportError:
		md5_new = None


class FastaFile(MutableMapping):

	def __init__(self,f,unmask=False,requireFai=False,ignoreFai=False,debug=None):
		"""
		f can be either a file object or a filename
		"""
		self.debug = debug if (debug != None) else []

		if (requireFai) and (ignoreFai):
			raise ValueError(".fai is required but is also to be ignored(!)")
		if (requireFai) and (type(f) != str):
			raise ValueError(".fai is required but fai filename cannot be constructed")

		fai = None
		if (type(f) == str) and ((requireFai) or (not ignoreFai)):
			if (f.endswith(".bgz")) or (f.endswith(".gz")) or (f.endswith(".bz2")) or (f.endswith(".zip")):
				raise ValueError("this .fai implementation doesn't support compressed fasta files (\"%d\")" % f)
			faiFilename = f+".fai"
			faiExists = exists(faiFilename)
			if (faiExists):
				fai = open(faiFilename,"r")
			elif (requireFai):
				raise ValueError(".fai is required but \"%d\" doesn't exist" % faiFilename)

		if (fai != None):
			f = open(f,"rb")
		elif (type(f) == str):
			f = open(f,"r")

		self.file   = f
		self.unmask = unmask

		self.seqCount = 0
		self.index = {}
		self.alias = {}

		if (fai == None):
			# read the sequences and create a name-to-sequence hash; we also
			# create an number-to-name hash; and if there is only one sequence,
			# we will map None to that name

			for (name,seq) in self.read_all_sequences(f):
				if (name in self.index):
					raise Exception("sequence name \"%s\" appears more than once in fasta file)" \
								  % (name,nameIx,nameOffset))
				if (unmask): seq = seq.upper()
				self.index[name] = seq
				self.alias[self.seqCount] = name
				self.seqCount += 1
				if ("signatures" in self.debug):
					seqHash = md5_new()
					seqHash.update(seq.encode("utf-8"))
					print("signature(\"%s\") = %s" % (name,seqHash.hexdigest().upper()),file=stderr)
			if (self.seqCount == 1): self.alias[None] = self.alias[0]
		else:
			# read the sequences and create a name-to-index-info hash (as well
			# as the number-to-name hash); entries in this hash will be
			# replaced by the sequence on an as-needed basis

			for (name,length,offset,lineBases,lineBytes) in self.read_fai(fai):
				if (name in self.index):
					raise Exception("sequence name \"%s\" appears more than once in fai file)" \
								  % (name,nameIx,nameOffset))
				self.index[name] = (length,offset,lineBases,lineBytes)
				self.alias[self.seqCount] = name
				self.seqCount += 1
			if (self.seqCount == 1): self.alias[None] = self.alias[0]

	def __getitem__(self,name):
		try:
			# convert to alias if necessary, and lookup sequence -or- fai info
			if (name not in self.index): name = self.alias[name]
			seq = self.index[name]
		except KeyError:
			return None
		if (type(seq) == tuple):
			# if we have fai info this sequence was not previously loaded; do
			# so now
			(length,offset,lineBases,lineBytes) = seq
			seq = self.read_sequence(self.file,length,offset,lineBases,lineBytes)
			if (self.unmask): seq = seq.upper()
			self.index[name] = seq
			if ("signatures" in self.debug):
				seqHash = md5_new()
				seqHash.update(seq.encode("utf-8"))
				print("signature(\"%s\") = %s" % (name,seqHash.hexdigest().upper()),file=stderr)
		return seq

	def __setitem__(self,name,val):
		raise Exception("FastaFile doesn't support __setitem__)")

	def __delitem__(self,name):
		raise Exception("FastaFile doesn't support __delitem__)")

	def __len__(self):
		return len(self.index)

	def __iter__(self):
		for name in self.index:
			yield name

	def keys(self):
		return self.index.keys()

	def read_all_sequences(self,f):
		seqName = None
		lineNum = 0
		for line in f:
			lineNum += 1
			line = line.strip()

			if (line.startswith(">")):
				if (seqName != None):
					seq = "".join(seq)
					yield (seqName,seq)
				fields = line[1:].split()
				if (fields == []):
					assert (False), \
					       "sequence has no name (at line %d)" % lineNum
				seqName = fields[0]
				seq = []
			elif (seqName == None):
				assert (False), "first sequence has no header"
			else:
				seq += [line]

		if (seqName != None):
			seq = "".join(seq)
			yield (seqName,seq)

	def read_sequence(self,f,length,offset,lineBases,lineBytes):
		seq = []
		lineOffset = offset
		basesLeft = length
		while (basesLeft > 0):
			f.seek(lineOffset)
			basesOnLine = lineBases if (basesLeft >= lineBases) else basesLeft
			seq += [f.read(basesOnLine).decode()]
			lineOffset += lineBytes
			basesLeft  -= basesOnLine
		return "".join(seq)

	def read_fai(self,f):
		# see http://www.htslib.org/doc/faidx.html
		lineNumber = 0
		for line in f:
			lineNumber += 1
			line = line.strip()

			fields = line.split()
			try:
				if (len(fields) != 5): raise ValueError
				name      = fields[0]
				length    = int(fields[1])
				offset    = int(fields[2])
				lineBases = int(fields[3])
				lineBytes = int(fields[4])
				if (length < 0): raise ValueError
				if (offset < 0): raise ValueError
				if (lineBases < 0): raise ValueError
				if (lineBytes < lineBases): raise ValueError
			except ValueError:
				raise ValueError("improper fai line (%d)\n%s" % (lineNumber,line))

			yield (name,length,offset,lineBases,lineBytes)

