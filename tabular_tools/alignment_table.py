#!/usr/bin/env python3
"""
Read pairwise alignments from lastz's general format. Typical input is below
but other columns may be included. It is up to the calling program what
columns are required.

    #name1    zstart1  end1     name2            strand2 zstart2+ end2+  id
    mule.chr1 9833952  9834016  donkey.contig12  -       142287   142351 97.7
    mule.chr1 8336423  8336491  donkey.contig12  +       30291    30359  95.9
    mule.chr1 25360782 25360844 donkey.contig12  -       160824   160886 98.6
    mule.chr1 11648921 11648978 donkey.contig163 -       51135    51192  95.9
    mule.chr1 11283997 11284065 donkey.contig163 -       57665    57733  99.0
    mule.chr1 18960023 18960091 donkey.contig163 +       56549    56617  96.9
     ...
"""

class Alignment: pass

class AlignmentTable(object):

	def __init__(self):
		"""
		We don't expet direct calls to this constructor. Instead, we expect
		objects to be created by calls to from_file().
		"""
		pass

	def __iter__(self):
		return self

	def __next__(self):
		if (self.iteratorStopped):
			if (self.filename != None): self.f.close()
			raise StopIteration
		elif (self.iteratorLookahead != None):
			a = self.iteratorLookahead
			self.iteratorLookahead = None
		else:
			try:
				a = self.iterator.__next__() # may raise StopIteration
			except StopIteration:
				if (self.filename != None): self.f.close()
				raise StopIteration
		return a

	@classmethod
	def from_file(cls,f,preFill=False,columnNames=None,
	              requiredColumns=None,nonRequiredColumns=None,columnAliases=None):
		"""
		Create an AlignmentTable from a file or stream.

		If preFill is True, we read the entire file into t.alignments.

		Otherwise (if preFill is False) we set t up as an iterator. The caller
		can then use the object like this:
		  t = AlignmentTable.from_file(...)
		  if (t.header != None):
		    do something with t.header
		  for a in t:
		    do something with alignment a

		f can either be a file object (already opened for read) or the name of
		a file (which we will open, read, and close).

		columnName, requiredColumns, nonRequiredColumns, and columnAliases are as
		defined in read_file().
		"""

		filename = None
		if (type(f) == str):
			filename = f
			if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
				from gzip import open as gzip_open
				f = gzip_open(filename,"rt")
			else:
				f = open(filename,"rt")

		t = AlignmentTable()
		t.f = f
		t.filename = filename
		t.header = None
		t.columnNames = columnNames

		# pre-fill the table

		if (preFill):
			t.iterator   = None
			t.alignments = []
			for a in AlignmentTable.read_file(f,columnNames=columnNames,
			                                  requiredColumns=requiredColumns,
			                                  nonRequiredColumns=nonRequiredColumns,
			                                  columnAliases=columnAliases):
				if (type(a) == str):
					t.header = a
				elif (type(a) == dict):
					if (t.columnNames == None):
						t.columnNames = a
				else:
					t.alignments += [a]
			if (filename != None):
				f.close()

		# don't pre-fill, set up an iterator

		else:
			t.alignments = None
			t.iterator   = AlignmentTable.read_file(f,columnNames=columnNames,
			                                        requiredColumns=requiredColumns,
			                                        nonRequiredColumns=nonRequiredColumns,
			                                        columnAliases=columnAliases)

			# fetch the first item, which may be the header

			try:
				a = t.iterator.__next__()
				t.iteratorStopped = False
			except StopIteration:
				t.iteratorStopped = True

			if (type(a) == str):
				t.header = a
				t.iteratorLookahead = None
			else:
				t.iteratorLookahead = a

			# fetch the second item, which may be the column names map

			if (not t.iteratorStopped):
				if (t.iteratorLookahead != None):
					a = t.iteratorLookahead
					t.iteratorLookahead = None
				else:
					try:
						a = t.iterator.__next__()
						t.iteratorStopped = False
					except StopIteration:
						t.iteratorStopped = True

				if (type(a) == dict):
					if (t.columnNames == None): t.columnNames = a
					t.iteratorLookahead = None
				else:
					t.iteratorLookahead = a

		return t


	@classmethod
	def read_file(cls,f,columnNames=None,requiredColumns=None,nonRequiredColumns=None,columnAliases=None):
		"""
		f is a file object (already opened for read).

		If columnNames is None, we require the first line of the file to be
		column headers, and extract the names from those. If columnNames is
		*not* None, the first line may or may not be column headers. If it is
		column headers, they must match what has been given in columnNames.
		columnNames is either a list of names, or a dict mapping a name to a
		column number (0 is the first column).

		requiredColumns lists the names of columns that are required to be
		present. The first six required columns *must* correspond to the fields
		name1, zstart1, end1, name2, zstart2+, and end2+. They don't have to
		have those specific names, but the must be conceptually equivalent to
		them.

		nonRequiredColumns lists the names of additional columns that, if
		present, will be carried as attributes in the alignment objects.

		columnAliases is a dict than maps an alias to an equivalent column
		name. For example, the input file might have a column "s" instead of
		"strand2".

		Note that we yield three types of objects:
		  - a string, which is the header line (if there is one)
		  - a dict, which maps names to column numbers (equivalent to
		    columnNames)
		  - a series of Alignments
		"""
		if (type(columnNames) == list):
			columnNames = {name:ix for (ix,name) in enumerate(columnNames)}

		if (requiredColumns == None):
			requiredColumns = ("name1","zstart1","end1",
			                   "name2","zstart2+","end2+",
			                   "strand")
		else:
			if (len(requiredColumns) < 6):
				raise ValueError("at least six columns must be specified in requiredColumns")
			requiredColumns = tuple(requiredColumns)

		if (columnAliases == None):
			columnAliases = {"strand2" : "strand",
			                 "s"       : "strand",
			                 "s2"      : "strand",
			                 "text1"   : "align1",
			                 "text2"   : "align2",
			                 "id%"     : "id",
			                 "con%"    : "con",
			                 "cov1%"   : "cov1",
			                 "cov2%"   : "cov2"}

		headerNames  = None
		name1Column  = None
		start1Column = None
		end1Column   = None
		name2Column  = None
		start2Column = None
		end2Column   = None

		numColumnsNeeded = None

		header = None
		lineNumber = 0
		for line in f:
			lineNumber += 1
			line = line.strip()

			# process the header line

			if (line.startswith("#")):
				if (header != None): continue  # we ignore all but the first header line
				header = line[1:]
				fields = line.split()
				fields[0] = fields[0][1:]
				headerNames = AlignmentTable.column_names(fields,requiredColumns,columnAliases=columnAliases)
				if (columnNames != None):
					differences = AlignmentTable.column_name_differences(headerNames,columnNames)
					if (differences != None):
						differences = "\n".join(["%s %d in file but spcified as %d"%(name,1+hIx,1+cIx) for (name,hIx,cIx) in differences])
						raise ValueError("column names provided within the file don't match those specified\n%s" % differences)
				else:
					columnNames = headerNames
					yield line
					yield columnNames
				continue

			# if we needed a header line and the input doesn't have one, it's
			# an error

			if (columnNames == None):
				raise ValueError("input column names weren't specified, and are not provided within the file")

			# on the arrival of the first non-header, derive some column
			# mapping info

			if (numColumnsNeeded == None):
				numColumnsNeeded = 1 + max([columnNames[name] for name in columnNames])

				(name1Column,start1Column,end1Column,
	             name2Column,start2Column,end2Column) = requiredColumns[:6]

				name1Column  = columnNames[name1Column ]
				start1Column = columnNames[start1Column]
				end1Column   = columnNames[end1Column  ]
				name2Column  = columnNames[name2Column ]
				start2Column = columnNames[start2Column]
				end2Column   = columnNames[end2Column  ]

			# split the line into fields

			fields = line.split()
			if (len(fields) < numColumnsNeeded):
				raise ValueError("not enough columns at line %d (%d, expected at least %d)" \
				               % (lineNumber,len(fields),numColumnsNeeded))

			# create the Alignment object and populate the common attributes

			a = Alignment()
			a.lineNumber = lineNumber
			a.line       = line
			a.name1      = fields[name1Column ]
			a.start1     = fields[start1Column]
			a.end1       = fields[end1Column  ]
			a.name2      = fields[name2Column ]
			a.start2     = fields[start2Column]
			a.end2       = fields[end2Column  ]

			# populate the remaining attributes

			for name in requiredColumns[6:]:
				if (name in columnAliases): name = columnAliases[name]
				val = fields[columnNames[name]]
				# $$$ perhaps recognize some fields and parse them, such as:
				# if (name == "score"): val = int_or_float(val)
				setattr(a,name,val)

			if (nonRequiredColumns != None):
				for name in nonRequiredColumns:
					if (name in columnAliases): name = columnAliases[name]
					if (name not in columnNames): continue
					val = fields[columnNames[name]]
					# $$$ perhaps recognize some fields and parse them, such as:
					# if (name == "score"): val = int_or_float(val)
					setattr(a,name,val)

			# validate the positional attributes
			# $$$ add a way to flip start2 and end2 for negative strand

			try:
				a.start1 = int(a.start1)
				a.end1   = int(a.end1)
				if (a.start1 >= a.end1): raise ValueError
			except ValueError:
				raise ValueError("bad alignment (at line %d), first species start/end\n%s" \
				               % (lineNumber,line))

			try:
				a.start2 = int(a.start2)
				a.end2   = int(a.end2)
				if (a.start2 >= a.end2): raise ValueError
			except ValueError:
				raise ValueError("bad alignment (at line %d), second species start/end\n%s" \
				               % (lineNumber,line))

			if ("strand" in requiredColumns):
				try:
					if (a.strand not in ["+","-"]): raise ValueError
				except ValueError:
					raise ValueError("bad alignment (at line %d), bad strand \"%s\"\n%s" \
					               % (lineNumber,a.strand,line))

			yield a

	@staticmethod
	def column_names(names,requiredColumns,columnAliases=None):
		if (columnAliases == None): columnAliases = {}
		columnNames = {}
		for (ix,name) in enumerate(names):
			actualName = name
			if (name in columnAliases): name = columnAliases[name]
			if (name not in requiredColumns): continue
			if (name in columnNames):
				raise ValueError(("\"%s\" (or alias) appears more than once in --format" % actualName))
			columnNames[name] = ix
		for name in requiredColumns:
			if (name not in columnNames):
				raise ValueError(("required name \"%s\" is absent" % name))
		return columnNames

	@staticmethod
	def column_name_differences(columnNames1,columnNames2):
		differences = []
		for name in columnNames2:
			if (name not in columnNames1):
				differences += [(name,-1,columnNames2[name])]
			if (columnNames1[name] != columnNames2[name]):
				differences += [(name,columnNames1[name],columnNames2[name])]
		for name in columnNames1:
			if (name not in columnNames2):
				differences += [(name,columnNames1[name],-1)]
		return None if (differences == []) else differences
