"""High-level PAM-scanning driver.

:func:`pamscan` is the single entry point shared by the command-line interface
(:mod:`pam_scanning.cli`) and the Tkinter GUI (:mod:`pam_scanning.gui`). It takes
the run parameters as keyword arguments, runs the full pipeline implemented in
:mod:`pam_scanning.library`, and writes a time-stamped results directory.
"""


def default_codon_table_path():
	"""Return the path to the yeast codon table bundled as package data."""
	from importlib import resources

	return str(resources.files("pam_scanning") / "data" / "codon_tables" / "yeast_64_1_1_all_nuclear.cusp.txt")


def pamscan(**kwargs):

	# Set **kwargs

	orf_file_path = kwargs['orf_file_path']
	orf_plus_buffer_file_path = kwargs['orf_plus_buffer_file_path']
	local_genome_file_path = kwargs['local_genome_file_path']
	codon_table_file_path = kwargs.get('codon_table_file_path', 'No file selected')
	codon_selection_file_path = kwargs['codon_selection_file_path']

	geneName = kwargs['geneName']
	localBlastDb = kwargs['localBlastDb']
	guidePrimerForwardSuffix = kwargs['guidePrimerForwardSuffix']
	insertPrimerForwardSuffix = kwargs['insertPrimerForwardSuffix']
	insertPrimerReverseSuffix = kwargs['insertPrimerReverseSuffix']

	primerLength = kwargs['primerLength']
	maxPamCutGap = kwargs['maxPamCutGap']
	codonsSamplingGap = kwargs['codonsSamplingGap']
	pamInclusionThreshold = kwargs['pamInclusionThreshold']
	pamInclusionSequenceThreshold = kwargs['pamInclusionSequenceThreshold'] 
	codonsSamplingGap = kwargs['codonsSamplingGap']

	outputPath = kwargs['outputPath']

	if orf_file_path == 'No file selected':
		print("ORF file is required")
		return 0

	if orf_plus_buffer_file_path == 'No file selected':
		print("ORF+ file is required")
		return 0

	if local_genome_file_path == 'No file selected':
		print("Local genome file is required")
		return 0

	# No codon table supplied: fall back to the bundled yeast codon table.
	if not codon_table_file_path or codon_table_file_path == 'No file selected':
		codon_table_file_path = default_codon_table_path()

	########################################################################################
	# Primer-order layout settings...
	########################################################################################

	# Microplate dimensions for the primer order...Options 96 or 384...
	microplateDimensions = 96

	########################################################################################
	# Get dependencies...
	########################################################################################

	from pam_scanning.library import fasta
	from pam_scanning.library import byRes, byCodon, singleLetterAA, reverseComplement
	from pam_scanning.library import findPamSites, tryToPamSilence, guideSilence, blastGuides, createPrimerOrder
	from pam_scanning.library import getOptiGuides, calculateScannableSequence
	from pam_scanning.library import set_codon_table
	import openpyxl

	set_codon_table(codon_table_file_path)

	# If they exist, collect user-provided codon selections...
	selectCodons = []
	if codon_selection_file_path and codon_selection_file_path != "No file selected":
		workbook = openpyxl.load_workbook(codon_selection_file_path)
		sheet = workbook["Mutations"]
		rows = sheet.max_row
		r = 2
		while r <= rows:
			cell = sheet.cell(row=r, column=1)
			residueNumber = cell.value
			# Necessary to skip previously occupied workbook cells that have been deleted... 
			if residueNumber == None:
				r += 1
				continue
			cell = sheet.cell(row=r, column=2)
			originalResidue = cell.value
			cell = sheet.cell(row=r, column=3)
			mutatedResidue = cell.value
			selectCodons.append((residueNumber, originalResidue, mutatedResidue))
			r += 1
	selectCodons.sort()

	# if not selectCodons:
	# 	print("No mutations were provided...abort PAM scan.")
	# 	raise SystemExit

	########################################################################################
	# Get DNA sequences...
	########################################################################################

	# A) Get the ORF sequence flanked by at least 100 bp of 5' and 3' genome sequence...
	sequenceInput = open(orf_plus_buffer_file_path, "r")
	orfPlusSequence = ""
	for line in sequenceInput.readlines():
		if line[0] == ">":
			continue
		if "\r\n" in line:
			sLine = line.split("\r\n") # Necessary becuase this is how SnapGene exports FASTA files...
		else:
			sLine = line.split("\n")
		if sLine:
			orfPlusSequence += sLine[0]
	sequenceInput.close()
	orfPlusSequence = orfPlusSequence.upper() # Initial sequence must be uppercase for functions to work properly...

	# B) Get the ORF sequence ATG to stop...
	sequenceInput = open(orf_file_path, "r")
	orfSequence = ""
	for line in sequenceInput.readlines():
		if line[0] == ">":
			continue
		if "\r\n" in line:
			sLine = line.split("\r\n") # Necessary becuase this is how SnapGene exports FASTA files...
		else:
			sLine = line.split("\n")
		if sLine:
			orfSequence += sLine[0]
	sequenceInput.close()
	orfSequence = orfSequence.upper() # Initial sequence must be uppercase for functions to work properly...

	########################################################################################
	# Make a time-stamped directory for the Pamscanning output....
	########################################################################################

	from os import mkdir, listdir, sep
	from os.path import exists
	from shutil import copyfile
	from time import strftime
	timeStamp = strftime("%Y.%m.%d-%H.%M.%S")

	# Root calculation results path...
	# outputPath = ".." + sep + "calculationResults" + sep
	outputPath += sep
	print("This is the output path", outputPath)
	if not exists(outputPath):
		mkdir(outputPath)
	outputPath += geneName + "-chimera-insertions-" + timeStamp + sep
	if not exists(outputPath):
		mkdir(outputPath)

	# QC path...
	qcPath = outputPath + "QC" + sep
	if not exists(qcPath):
		mkdir(qcPath)

	# BLAST+ path...
	blastPath = outputPath + "BLAST+" + sep
	if not exists(blastPath):
		mkdir(blastPath)

	# ORDER path...
	orderPath = outputPath + "ORDER" + sep
	if not exists(orderPath):
		mkdir(orderPath)

	if not exists(qcPath + "solutionsFasta" + sep):
		mkdir(qcPath + "solutionsFasta" + sep)
	copyfile(orf_file_path, qcPath + orf_file_path.split(sep)[-1])

	########################################################################################
	# Copy the orf plus sequence to file for quality control checks...
	########################################################################################

	copyfile(orf_plus_buffer_file_path, qcPath + orf_plus_buffer_file_path.split(sep)[-1])

	########################################################################################
	# Set codon information...
	########################################################################################

	orfStartIndex = orfPlusSequence.index(orfSequence)
	orfEndIndex = orfPlusSequence.index(orfSequence) + len(orfSequence)
	orfOffset = orfPlusSequence.index(orfSequence)
	i, allCodons = orfPlusSequence.index(orfSequence), {}
	while i < orfEndIndex:
		codon = orfPlusSequence[i:i+3]
		allCodons[(i, i+1, i+2)] = (codon, int(((i-orfOffset)/3) + 1))
		i += 3

	# Filter codons by A) selectCodons or B) codonSamplingGap if set...
	codonKeys = sorted(allCodons)
	# A) Filter by selectCodons...
	filteredCodons = {}
	if selectCodons:
		for selectCodon in selectCodons:
			residueNumber = selectCodon[0]
			for codonKey in codonKeys:
				if residueNumber == allCodons[codonKey][1]:
					filteredCodons[codonKey] = allCodons[codonKey]
					break
	# B) Filter by codonSamplingGap
	else:
		i = 0
		while i < len(codonKeys):
			codonKey = codonKeys[i]
			filteredCodons[codonKey] = allCodons[codonKey]
			i += codonsSamplingGap
	codons = filteredCodons
	codonKeys = sorted(codons)

	########################################################################################
	# Identify and attempt to silence PAM sites...
	########################################################################################

	result = findPamSites(orfPlusSequence, orfStartIndex, orfStartIndex + len(orfSequence))
	guidesF, guidesR = result[0], result[1]
	guides = guidesF
	guides.update(guidesR)
	result = tryToPamSilence(orfPlusSequence, orfStartIndex, orfEndIndex, guides)
	silencedGuides, unsilencedGuides = result[0], result[1]

	########################################################################################
	# Attempt to GUIDE silence any guide that cannot be PAM silenced...
	########################################################################################

	guideSilencers = guideSilence(orfPlusSequence, orfStartIndex, unsilencedGuides)
	silencedGuides.update(guideSilencers)
	for key in silencedGuides:
		if key in unsilencedGuides:
			del unsilencedGuides[key]
	# Update guides to only silenceable guides...
	tmpGuides = {}
	keys = sorted(silencedGuides)
	for key in keys:
		tmpGuides[key] = guides[key]
	guides = tmpGuides

	########################################################################################
	# BLAST guides against the reference genome to avoid off-target Cas9 cutting...
	########################################################################################

	result = blastGuides(blastPath, geneName, guides, localBlastDb, local_genome_file_path, pamInclusionSequenceThreshold=pamInclusionSequenceThreshold, pamInclusionThreshold=pamInclusionThreshold)
	safeGuides, unsafeGuides, pamInclusionsDict, superConservativePamInclusionDict = result[0], result[1], result[2], result[3]
	keys, tmpGuides, tmpSilencedGuides = sorted(safeGuides), {}, {}
	for key in keys:
		# Remove any guides that could lead to off-target cutting...
		tmpGuides[key] = guides[key]
		tmpSilencedGuides[key] = silencedGuides[key]
	guides, silencedGuides = tmpGuides, tmpSilencedGuides

	# Identify the optimal safeguide for each codon (post BLAST verification)....
	# This information is used to calculate total PAM-scannable sequence for the ORF...
	optiGuidesAllCodons = getOptiGuides(allCodons, guides, orfPlusSequence, maxPamCutGap, pamInclusionsDict)
	keys, tmpGuides, tmpSilencedGuides = sorted(optiGuidesAllCodons), {}, {}
	for key in keys:
		if not optiGuidesAllCodons[key]:
			continue
		guideKey = optiGuidesAllCodons[key][0]
		guideSeq = optiGuidesAllCodons[key][1]
		guideCut = optiGuidesAllCodons[key][2]
		tmpGuides[guideKey] = guides[guideKey]
		tmpSilencedGuides[guideKey] = silencedGuides[guideKey]
	guides, silencedGuides = tmpGuides, tmpSilencedGuides

	# Identify the optimal safeguide for each codon of interest (post BLAST verification)....
	optiGuides = getOptiGuides(codons, guides, orfPlusSequence, maxPamCutGap, pamInclusionsDict)
	keys, noOptiGuideCodonSet = sorted(optiGuides), {}
	for key in keys:
		if not optiGuides[key]:
			noOptiGuideCodonSet[key] = None

	########################################################################################
	# Make and bundle necessary guide and payload primers...
	########################################################################################

	wb = openpyxl.Workbook()
	sheet = wb["Sheet"]
	wb.remove(sheet)
	sheet = wb.create_sheet("Guide Solutions By Row")
	cell = sheet.cell(row=1, column=1)
	cell.value = "Codon"
	cell = sheet.cell(row=1, column=2)
	cell.value = "Cut @ Base"
	cell = sheet.cell(row=1, column=3)
	cell.value = "PAM Cut Gap"
	cell = sheet.cell(row=1, column=4)
	cell.value = "PAM Inclusions"
	cell = sheet.cell(row=1, column=5)
	cell.value = "Original Guide"
	cell = sheet.cell(row=1, column=6)
	cell.value = "Silenced Guide"
	cell = sheet.cell(row=1, column=7)
	cell.value = "Primer"
	cell = sheet.cell(row=1, column=8)
	cell.value = "Primer Name"
	cell = sheet.cell(row=1, column=9)
	cell.value = "Primer Guide"

	codonKeys, row, primerOrderGuides, primerOrderInserts = sorted(optiGuides), 2, {}, {}
	for codonKey in codonKeys:

		if codonKey in noOptiGuideCodonSet:
			continue

		# Get codon information...
		codon = codons[codonKey]
		codonBases = codon[0]
		codonNumber = codon[1]

		# Build primers...
		guideKey = optiGuides[codonKey][0]
		guideSeq = optiGuides[codonKey][1].lower()
		guideGap = optiGuides[codonKey][2]
		guideInclusions = optiGuides[codonKey][3]
		silencedGuide = silencedGuides[guideKey][0]
		silencedOrf = silencedGuides[guideKey][1]
		primerLengthF = primerLength - len(insertPrimerForwardSuffix)
		primerLengthR = primerLength - len(insertPrimerReverseSuffix)
		leftHomologyArm = silencedOrf[codonKey[0]-primerLengthF+3:codonKey[0]+3]
		rightHomologyArm = silencedOrf[codonKey[0]+3:codonKey[0]+3+primerLengthR]
		primerF = leftHomologyArm + insertPrimerForwardSuffix
		primerR = reverseComplement(rightHomologyArm) + insertPrimerReverseSuffix
		primerNameGuide = str(guideKey[0]) + ".gF"
		primerNameInsertionF = str(guideKey[0]) + ".iF." + str(codonNumber)
		primerNameInsertionR = str(guideKey[0]) + ".iR." + str(codonNumber)
		if guideSeq[-2:] == "gg":
			primerOrderGuides[primerNameGuide] = guideSeq[:-3] + guidePrimerForwardSuffix
		else:
			primerOrderGuides[primerNameGuide] = reverseComplement(guideSeq)[:-3] + guidePrimerForwardSuffix
		primerOrderInserts[primerNameInsertionF] = primerF
		primerOrderInserts[primerNameInsertionR] = primerR

		# Write .fa file for viewing in SnapGene...
		fastaHeader = "> "
		fastaHeader += "Codon: " + str(codonNumber)
		fastaHeader += " | Original guide: " + guideSeq 
		fastaHeader += " | Silenced guide: " + silencedGuide
		fastaHeader += " | Left homology arm: " + leftHomologyArm
		fastaHeader += " | Right homology arm: " + rightHomologyArm 
		fastaHeader += " | Gap: " + str(guideGap)  
		fastaHeader += " | Inclusions: " + str(guideInclusions) + "\n"
		output = open(qcPath + "solutionsFasta" + sep + geneName + "-" + str(codonNumber) + ".fa", "w")
		output.write(fastaHeader)
		output.write(fasta(silencedOrf))
		output.close()

		# Add solution to the .xls solution workbook...
		# Forward...
		cell = sheet.cell(row=row, column=1)
		cell.value = codonNumber
		cell = sheet.cell(row=row, column=2)
		cell.value = guideKey[0]
		cell = sheet.cell(row=row, column=3)
		cell.value = guideGap
		cell = sheet.cell(row=row, column=4)
		cell.value = guideInclusions
		cell = sheet.cell(row=row, column=5)
		cell.value = guideSeq
		cell = sheet.cell(row=row, column=6)
		cell.value = silencedGuide
		cell = sheet.cell(row=row, column=7)
		cell.value = primerF
		cell = sheet.cell(row=row, column=8)
		cell.value = primerNameInsertionF
		cell = sheet.cell(row=row, column=9)
		cell.value = primerNameGuide

		# Reverse...
		cell = sheet.cell(row=row+1, column=1)
		cell.value = codonNumber 
		cell = sheet.cell(row=row+1, column=2)
		cell.value = guideKey[0]
		cell = sheet.cell(row=row+1, column=3)
		cell.value = guideGap
		cell = sheet.cell(row=row+1, column=4)
		cell.value = guideInclusions
		cell = sheet.cell(row=row+1, column=5)
		cell.value = guideSeq
		cell = sheet.cell(row=row+1, column=6)
		cell.value = silencedGuide
		cell = sheet.cell(row=row+1, column=7)
		cell.value = primerR
		cell = sheet.cell(row=row+1, column=8)
		cell.value = primerNameInsertionR
		cell = sheet.cell(row=row+1, column=9)
		cell.value = primerNameGuide

		row += 2

	wb.save(qcPath + geneName + "-guideSolutions.xlsx")

	########################################################################################
	# BLAST guides against the reference genome...
	# At the this point, the remaining guides are safe and avoid off-target Cas9 cutting...
	# The purpose of this code is to generate the x-guidesBlastResult-final.txt
	########################################################################################

	optiGuidesForFinalBlast = {}
	keys = sorted(optiGuides)
	for key in keys:

		if not optiGuides[key]:
			continue

		guideKey = optiGuides[key][0]
		guideSequence = guideKey[1]
		optiGuidesForFinalBlast[guideKey] = guideSequence

	result = blastGuides(blastPath, geneName, optiGuidesForFinalBlast, localBlastDb, local_genome_file_path, pamInclusionSequenceThreshold=pamInclusionSequenceThreshold, pamInclusionThreshold=pamInclusionThreshold, final=1)
	pamInclusionsDict, superConservativePamInclusionDict = result[2], result[3]

	########################################################################################
	# Create the primer order...
	########################################################################################

	createPrimerOrder(orderPath, microplateDimensions, geneName, primerOrderGuides, primerOrderInserts)

	########################################################################################
	# Report PAM inclusions warning for potential off-target cutting...
	########################################################################################

	codonKeys, writeInclusions, writeInclusionsSuper = sorted(optiGuides), {}, {}
	for codonKey in codonKeys:

		if codonKey in noOptiGuideCodonSet:
			continue

		# Get codon information...
		codon = codons[codonKey]
		codonBases = codon[0]
		codonNumber = codon[1]
		optiGuide = optiGuides[codonKey][0][1]

		# Conservative PAM inclusions...
		keys = sorted(pamInclusionsDict)
		for key in keys:
			inclusions = pamInclusionsDict[key]
			inclusionGuide = inclusions[0]
			if optiGuide == inclusionGuide:
				print("Warning: Some guides will contain potential PAM inclusions...", inclusions)
				writeInclusions[(codonNumber, key)] = pamInclusionsDict[key]

		# Super conservative PAM inclusions...
		keys = sorted(superConservativePamInclusionDict)
		for key in keys:
			inclusions = superConservativePamInclusionDict[key]
			inclusionGuide = inclusions[0]
			if optiGuide == inclusionGuide:
				writeInclusionsSuper[(codonNumber, key)] = superConservativePamInclusionDict[key]

	if writeInclusions:

		# Warning path...
		warningPath = outputPath + "WARNINGS/"
		if not exists(warningPath):
			mkdir(warningPath)

		f = open(warningPath + geneName + "-pamInclusionWarnings.txt", "w")
		sFormat = "%-6s | %-25s | %-25s | %-50s\n"
		f.write(sFormat % ("Codon", "Guide", "BLAST hit", "BLAST result"))
		sFormat = "%-6i | %-25s | %-25s | %-50s\n"
		for key in sorted(writeInclusions):
			inclusions = writeInclusions[key]
			f.write(sFormat % (key[0], inclusions[0], inclusions[1], inclusions[2]))
		f.close()

	# Don't write the super conservative inclusion file if it is the same as the parameter-defined inclusion file...
	if writeInclusionsSuper and len(superConservativePamInclusionDict) != len(pamInclusionsDict):

		# Warning path...
		warningPath = outputPath + "WARNINGS/"
		if not exists(warningPath):
			mkdir(warningPath)

		f = open(warningPath + geneName + "-pamInclusionWarnings-superConservative.txt", "w")
		sFormat = "%-6s | %-25s | %-25s | %-50s\n"
		f.write(sFormat % ("Codon", "Guide", "BLAST hit", "BLAST result"))
		sFormat = "%-6i | %-25s | %-25s | %-50s\n"
		for key in sorted(writeInclusionsSuper):
			inclusions = writeInclusionsSuper[key]
			f.write(sFormat % (key[0], inclusions[0], inclusions[1], inclusions[2]))
		f.close()

	if noOptiGuideCodonSet:

		# Warning path...
		warningPath = outputPath + "WARNINGS/"
		if not exists(warningPath):
			mkdir(warningPath)

		f = open(warningPath + geneName + "-unscannableCodons.txt", "w")
		keys = sorted(noOptiGuideCodonSet)
		for key in keys:
			# Get codon information...
			codon = codons[key]
			codonBases = codon[0]
			codonNumber = codon[1]
			f.write(str(codonNumber) + "\n")
		f.close()

	########################################################################################
	# Report the total PAM-scannable sequence of the ORF for reference...
	########################################################################################

	calculateScannableSequence(qcPath, geneName, orfSequence, orfPlusSequence, safeGuides, maxPamCutGap)


