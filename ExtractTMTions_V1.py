#!/usr/bin/python3
# TMTpro Intensity Extractor and Plink2 Mapping

#used to get cmd line arguments 
import argparse, os, glob, re, math, statistics


#Read command line arguments and create help documentation using argparse
parser = argparse.ArgumentParser(
	description='''Tool for extracting TMT ion intensities and optionally mapping to plink2 output. Dr James Wright - 2024 - Institute of Cancer Research''')

parser.add_argument('--mgfdir', '-i', dest='mgfdir', help='Directory conatining MGF spectrum files matching the plink2 input', required=True)
parser.add_argument('--tol', '-t', dest='itol', default='15',  help='TMT ion matching tolerance inppm, default = 15ppm')
parser.add_argument('--proteins', '-p', dest='proteinfile', default='', help='csv plink2 output File Containing proteins and peptides')
parser.add_argument('--crosslinks', '-c', dest='crosslinkfile', default='', help='csv plink2 output File Containing crosslinked peptides')
parser.add_argument('--outprefix', '-o', dest='outpre', default='', help='prefix for output files')
parser.add_argument('--norm', '-n', dest='normal', default=False, action='store_true',  help='Run normalisation and scaling of TMT channels using total TMT channel intensity (beta method). Default=False ')

args = parser.parse_args()


#Dictionary of TMT ions to extract
TMTproIons = {
	'126' : 126.127726,
	'127N' : 127.124761,
	'127C' : 127.131081,
	'128N' : 128.128116,
	'128C' : 128.134436,
	'129N' : 129.131471,
	'129C' : 129.13779,
	'130N' : 130.134825,
	'130C' : 130.141145,
	'131N' : 131.13818,
	'131C' : 131.1445,
	'132N' : 132.141535,
	'132C' : 132.147855,
	'133N' : 133.14489,
	'133C' : 133.15121,
	'134N' : 134.148245,
	'134C' : 134.154565,
	'135N' : 135.1516
}


######### Function ExtractTMTfromMGF ############
#
# This function parses MGF formatted spectra and extracts TMT ion intensities that match with in tolerence x
#
#
def MGFExtractTMTions(mgf, ions, ppm, TMTpeaks, norm):

	fid = os.path.basename(mgf)
	sid = ''

	TMTsum = {}
	TMTpeaks[fid] = {}

	#Covert PPM to Da for each TMT ion to match
	#PPM = ( theoMZ - expMZ / theoMZ) * 1000000
	tols = {}
	for i in ions:

		TMTsum[i] = 0

		tols[i] = {}
		da = (int(ppm) / 1000000) * ions[i]

		tols[i]['min'] = ions[i] - da
		tols[i]['max'] = ions[i] + da

	MGF = open(mgf, 'r')

	#loop each line in the file
	for line in MGF:
		#Extract Spectrum Title
		if line.startswith('TITLE='):
			sid = line.rstrip().replace('TITLE=', '')
			TMTpeaks[fid][sid] = {}

		#If spectral peak line
		if line[0].isdigit():
			ln = line.rstrip().split()
			mz = float(ln[0])

			#If in TMT mz range
			if 125 < mz < 136:

				for i in ions:
					if tols[i]['min'] <= mz <= tols[i]['max']:

						err =  int ( ( ( ions[i] - mz ) / ions[i] ) * 1000000 )

						TMTsum[i] += float(ln[1])

						if i not in TMTpeaks[fid][sid]:
							TMTpeaks[fid][sid][i] = {}
							TMTpeaks[fid][sid][i]['int'] = float(ln[1])
							TMTpeaks[fid][sid][i]['err'] = err
						else :
							print ("Warning: TMT Peak " + str(i) + " is not unique in spectrum " + str(sid) + " in MGF " + str(fid) + " at tolerance " + str(ppm) + " with error (ppm) " + str(err) + " !!!")
							TMTpeaks[fid][sid][i]['int'] += float(ln[1])
							if abs(err) > abs(TMTpeaks[fid][sid][i]['err']):
								TMTpeaks[fid][sid][i]['err'] = err

	MGF.close()

	
	if norm:
		#NORMALISE RAW INTENSITIES
		maxchannel = max(TMTsum, key=TMTsum.get)
		for i in ions:
			normfactor = TMTsum[maxchannel] / TMTsum[i]
			print ("Channel: " + str(i) + " has normalisation factor " + str(normfactor))
			if normfactor >= 100:
				print ("Warning: Channel " + str(i) + " has normalisation factor greater than 100 (" + str(normfactor) + ") and is being excluded from normalisation and scaling")
			else :
				for sid in TMTpeaks[fid]:
					if i in TMTpeaks[fid][sid]:
						TMTpeaks[fid][sid][i]['norm'] = TMTpeaks[fid][sid][i]['int'] * normfactor

		#SCALE INTENSITIES
		for sid in TMTpeaks[fid]:
			summed = 0
			count = 0
			for i in ions:
				if i in TMTpeaks[fid][sid]:
					if 'norm' in TMTpeaks[fid][sid][i]:
						summed += TMTpeaks[fid][sid][i]['norm']
						count += 1

			if count > 0 and summed > 0:
				mean = summed / count
				for i in ions:
					if i in TMTpeaks[fid][sid]:
						if 'norm' in TMTpeaks[fid][sid][i]:
							TMTpeaks[fid][sid][i]['scaled'] = (TMTpeaks[fid][sid][i]['norm'] / mean) * 100

	return


##################### END FUNCTION ############################

######### Function GenerateTMTionTable ############
#
# Writes a file containing all the extracted TMT ions and intensities matched in MGF
#
def WriteTMTionTable(out, ions, TMTpeaks, norm):

	outFile =out + "_TMTionTable.txt"
	outTMT = open(outFile, 'w')

	outTMT.write("MGF\tTITLE")
	for i in sorted(ions):
		outTMT.write("\t" + i + "-intensity_raw")
	for i in sorted(ions):
		outTMT.write("\t" + i + "-err")

	if norm:
		for i in sorted(ions):
			outTMT.write("\t" + i + "-intensity_normalised")
		for i in sorted(ions):
			outTMT.write("\t" + i + "-intensity_scaled")

	outTMT.write("\n")

	for f in TMTpeaks:
		for s in TMTpeaks[f]:
			if 'MEDIAN' not in s:
				outTMT.write( str(f) + "\t" + str(s) )
				for i in sorted(ions):
					if i not in TMTpeaks[f][s]:
						outTMT.write( "\t-" )
					else:
						outTMT.write( "\t" + str(TMTpeaks[f][s][i]['int']) )

				for i in sorted(ions):
					if i not in TMTpeaks[f][s]:
						outTMT.write( "\t-" )
					else:
						outTMT.write( "\t" + str(TMTpeaks[f][s][i]['err']) )

				if norm:
					for i in sorted(ions):
						if i not in TMTpeaks[f][s]:
							outTMT.write( "\t-" )
						else:
							if 'norm' in TMTpeaks[f][s][i]:
								outTMT.write( "\t" + str(round(TMTpeaks[f][s][i]['norm'], 2)) )
							else:
								outTMT.write( "\t-" )

					for i in sorted(ions):
						if i not in TMTpeaks[f][s]:
							outTMT.write( "\t-" )
						else:
							if 'scaled' in TMTpeaks[f][s][i]:
								outTMT.write( "\t" + str(round(TMTpeaks[f][s][i]['scaled'], 1)) )
							else:
								outTMT.write( "\t-" )
						
				outTMT.write("\n")

	outTMT.close()

	return

##################### END FUNCTION ############################

######### Function ReWrite Plink Results ############
#
# Parses a PLink output file and annotates TMT ion intensities extracted from MGF into the file
# 
#
def MapTMTionsToPlink(out, file, ions, TMTpeaks, norm):

	parent_line = ""
	psm_lines = ""
	parent_stats = {}
	psm_count = 0

	outFile=out + "_TMTmapped_" + os.path.basename(file)
	with open(outFile, 'w') as outPLINK:

		INP = open(file, 'r')
		for line in INP:

			#PROTEIN/PEPTIDE HEADER
			if re.match("[^,]+_Order", line): 
				outPLINK.write( line.rstrip() + ",PSM_Count" )
				for i in sorted(ions):
					outPLINK.write("," + i + "-Sum_intensity")

				if norm:
					for i in sorted(ions):
						outPLINK.write("," + i + "-Sum_intensity_normalised")
					for i in sorted(ions):
						outPLINK.write("," + i + "-Sum_intensity_scaled")

				outPLINK.write("\n")

			#PSM HEADER
			elif re.match(",Spectrum_Order", line): 
				outPLINK.write( line.rstrip() )
				outPLINK.write(",MGF,TITLE")
				for i in sorted(ions):
					outPLINK.write("," + i + "-intensity")

				if norm:
					for i in sorted(ions):
						outPLINK.write("," + i + "-intensity_normalised")
					for i in sorted(ions):
						outPLINK.write("," + i + "-intensity_scaled")

				outPLINK.write("\n")

			#PROTEIN/PEPTIDE DATA
			elif re.match("\d+", line): 

				if parent_line:
					outPLINK.write( parent_line + "," + str(psm_count))
					## CALCULATE AND ADD PROTEIN/PEPTIDE STATS
					if psm_count > 0:

						nSum = 0
						nCount = 0
						for i in sorted(ions):
							outPLINK.write( "," + str(parent_stats[i]['int']))

						if norm:
							for i in sorted(ions):
								outPLINK.write( "," + str(round(parent_stats[i]['norm'],2)))
								if parent_stats[i]['norm'] > 0:
									nSum += parent_stats[i]['norm']
									nCount += 1

							if nCount > 0:
								nMean = nSum / nCount
								for i in sorted(ions):
									if parent_stats[i]['norm'] > 0:
										scaled = (parent_stats[i]['norm'] / nMean) * 100
										outPLINK.write( "," + str(round(scaled,1)))
									else:
										outPLINK.write( ",-")
							else:
								for i in sorted(ions):
									outPLINK.write( ",-")

					else:
						for i in sorted(ions):
							outPLINK.write( ",-")

						if norm:
							for i in sorted(ions):
								outPLINK.write( ",-")
							for i in sorted(ions):
								outPLINK.write( ",-")

					outPLINK.write( "\n" )
					outPLINK.write( psm_lines )

				#Reset peptide/protein/psm 
				parent_line = line.rstrip() 
				psm_lines = ""
				parent_stats = {}
				psm_count = 0

				for i in ions:
					parent_stats[i] = {}
					parent_stats[i]['int'] = 0
					parent_stats[i]['norm'] = 0

			#PSM DATA
			elif re.match(",\d+", line): 
				psm_lines +=( line.rstrip() )
				cols = line.rstrip().split(',')	
				match = False
				psmquantified = False
				for f in TMTpeaks:
					if cols[2] in TMTpeaks[f]:
						if match:
							print ("Warning: Multiple Spectrum Match for PSM in file " + str(f) + "!")
							print (line)

						s = cols[2]
						match = True
						psm_lines +=( "," + str(f) + "," + str(s) )
						for i in sorted(ions):
							if i not in TMTpeaks[f][s]:
								psm_lines +=( ",-" )
							else:
								psm_lines +=( "," + str(TMTpeaks[f][s][i]['int']) )
								parent_stats[i]['int'] += TMTpeaks[f][s][i]['int']

						if norm:
							for i in sorted(ions):
								if i not in TMTpeaks[f][s]:
									psm_lines +=( ",-" )
								else:
									if 'norm' in TMTpeaks[f][s][i]:
										psm_lines +=( "," + str(round(TMTpeaks[f][s][i]['norm'], 2)) )
										parent_stats[i]['norm'] += TMTpeaks[f][s][i]['norm']
										psmquantified = True
									else:
										psm_lines +=( ",-" )

							for i in sorted(ions):
								if i not in TMTpeaks[f][s]:
									psm_lines +=( ",-" )
								else:
									if 'scaled' in TMTpeaks[f][s][i]:
										psm_lines +=( "," + str(round(TMTpeaks[f][s][i]['scaled'], 1)) )
									else:
										psm_lines +=( ",-" )

				if match is False:
					print ("Warning: No Spectrum Match for PSM!")
					print (line)
					psm_lines +=( ",-,-")
					for i in sorted(ions):
						psm_lines +=( ",-")
					if norm:
						for i in sorted(ions):
							psm_lines +=( ",-")
						for i in sorted(ions):
							psm_lines +=( ",-")

				psm_lines +=("\n")

				if psmquantified:
					psm_count += 1

			else:
				outPLINK.write(line)	

		if parent_line:
			outPLINK.write( parent_line + "," + str(psm_count))
			## CALCULATE AND ADD PROTEIN/PEPTIDE STATS
			if psm_count > 0:

				nSum = 0
				nCount = 0
				for i in sorted(ions):
					outPLINK.write( "," + str(parent_stats[i]['int']))

				if norm:
					for i in sorted(ions):
						outPLINK.write( "," + str(round(parent_stats[i]['norm'],2)))
						if parent_stats[i]['norm'] > 0:
							nSum += parent_stats[i]['norm']
							nCount += 1

					if nCount > 0:
						nMean = nSum / nCount
						for i in sorted(ions):
							if parent_stats[i]['norm'] > 0:
								scaled = (parent_stats[i]['norm'] / nMean) * 100
								outPLINK.write( "," + str(round(scaled,1)))
							else:
								outPLINK.write( ",-")
					else:
						for i in sorted(ions):
							outPLINK.write( ",-")

			else:
				for i in sorted(ions):
					outPLINK.write( ",-")

				if norm:	
					for i in sorted(ions):
						outPLINK.write( ",-")
					for i in sorted(ions):
						outPLINK.write( ",-")

			outPLINK.write( "\n" )
			outPLINK.write( psm_lines )	
	
	outPLINK.close()

	return

##################### END FUNCTION ############################

######## MAIN ##########

print ("TMT Ion Extractor and Plink2 Mapper V1")
print ("Dr James Wright - ICR - 2023")


#1. Extract TMT Ion Intensities from each MGF file
print ("1. Starting Extraction of ions from MGF...")
TMTpeaks = {}
#Loop through each MGF file in results directory
for mgf in glob.glob(args.mgfdir + '/*.mgf'):
	print ("Found MGF: " + mgf)
	MGFExtractTMTions(mgf, TMTproIons, args.itol, TMTpeaks, args.normal)
print ("\t..done")


#2. Write TMT Ion Intensities to file
print ("2. Writing ion table...")
WriteTMTionTable(args.outpre, TMTproIons, TMTpeaks, args.normal)
print ("\t..done")


#3. Map to Plink Protein File
if args.proteinfile:
	print ("3. Mapping Ions to Protein Results...")
	MapTMTionsToPlink(args.outpre, args.proteinfile, TMTproIons, TMTpeaks, args.normal)
	print ("\t..done")

#4. Map to Plink CrossLinker File
if args.crosslinkfile:
	print ("4. Mapping Ions to Protein Results...")
	MapTMTionsToPlink(args.outpre, args.crosslinkfile, TMTproIons, TMTpeaks, args.normal)
	print ("\t..done")

print ("Analysis Complete!")


#######################

