#!/usr/bin/env python
# encoding: utf-8

import re
import sys
import gzip
import pysam
#import itertools
import mimetypes
from collections import OrderedDict  #,defaultdict


class VCF(object):
    """docstring for VCF"""

    def __init__(self,
                 input,
                 output=None,
                 populations=None,
                 region=None,
                 window_size=1,
                 step=0,
                 snvs=None):

        super(VCF, self).__init__()

        self.input = input
        self.output = output
        self.step = step
        self.populations = populations
        self.region = region
        self.window_size = window_size
        self.header = self.__extract_header__()
        self.chrms_2_sizes = self.__get_chrm_ids_and_sizes__()
        self.empty_vcf_line = self.make_empty_vcf_ordered_dict()


    def __open_vcf__(self):
        """Open vcf file as gzip or as text."""

        if type(self.input) is file:
            fin = self.input

        elif mimetypes.guess_type(self.input)[-1] is 'gzip':
            fin = gzip.open(self.input, 'rb')

        else:
            fin = open(self.input, 'r')

        return fin

    def __extract_header__(self):
        #print [line.strip() for line in self.__open_vcf__()]

        header = []
        for line in self.__open_vcf__():

            if line.startswith("#") == True:
                header.append(line.strip())
            else:
                break

        return header


    def __get_chrm_ids_and_sizes__(self):
        """ Extract chromosome ids and sizes from vcf file.
            Return as dictionary"""

        chrms_sizes_dict = OrderedDict()
        # with self.__open_vcf__() as fin:

        for line in self.header:

            if line.startswith("##contig"):
                chrm_name = re.findall(r'ID=.*,', line)
                chrm_name = chrm_name[0].strip('ID=').strip(',')

                chrm_length = re.findall(r'length=.*>', line)
                chrm_length = int(chrm_length[0].strip('length=').strip('>'))

                chrms_sizes_dict[chrm_name] = chrm_length
                #break

        return chrms_sizes_dict

    def make_empty_vcf_ordered_dict(self):
        """Open VCF file and read in #CHROM line as an Ordered Dict"""

        header_dict = None
        for line in self.header:
            if line.startswith("#CHROM"):
                header = line.strip("#").strip().split()
                header_dict = OrderedDict([(item, None) for item in header])
                break

        return header_dict


    def process_snp(self, snp_call):
        if snp_call == "0/0":
            return (0,0)

        elif snp_call == "1/1":
            return (1,1)

        elif snp_call == '1/0' or \
               snp_call == '0/1':
            return (0,1)

        # skip multiallelic sites
        else:
            return None


    def process_snp_call(self, snp_call, ref, alt, IUPAC_ambiguities=False):
        """Process VCF genotype fields."""

        # IUPAC ambiguity codes
        IUPAC_dict = {('A', 'C'): 'M',
                      ('A', 'G'): 'R',
                      ('A', 'T'): 'W',
                      ('C', 'G'): 'S',
                      ('C', 'T'): 'Y',
                      ('G', 'T'): 'K',
                      ('A', 'C', 'G'): 'V',
                      ('A', 'C', 'T'): 'H',
                      ('A', 'G', 'T'): 'D',
                      ('C', 'G', 'T'): 'B'}


        if snp_call == None:
            called_base = "-"

        if snp_call["GT"] == "./.":
            called_base = "-"

        else:
            allele1, allele2 = snp_call["GT"].split("/")

            # process "0/0"
            if allele1 == '0' and allele2 == '0':
                called_base = ref

            if allele1 == '1' and allele2 == '1':
                called_base = alt

            # process "0/N"
            if allele1 == '0' and allele2 != '0':

                if IUPAC_ambiguities == False:
                    called_base = 'N'

                else:
                    call = [ref] + [alt.split(',')[int(allele2) - 1]]
                    call.sort()
                    call = tuple(call)
                    called_base = IUPAC_dict[call]

            # process "2/2, 1/2, etc."
            if int(allele1) >= 1 and int(allele2) > 1:

                # deal with homozygotes
                if allele1 == allele2:
                    called_base = alt.split(',')[int(allele1) - 1]

                # deal with heterozygotes
                else:

                    if IUPAC_ambiguities == False:
                        called_base = 'N'

                    else:
                        ref = alt.split(',')[int(allele1) - 1]
                        alt = alt.split(',')[int(allele2) - 1]
                        call = [ref, alt]
                        call.sort()
                        call = tuple(call)
                        called_base = IUPAC_dict[call]

        return called_base

    def count_alleles(self, chunk):

        results = []
        for line in chunk:
            pop_counts = {}

            for pop in self.populations.keys():
                allele_counts = {'REF':0, 'ALT':0}

                for sample in self.populations[pop]:

                    if line[sample] != None:
                        ref, alt = self.process_snp(line[sample]['GT'])
                        allele_counts['REF'] += ref
                        allele_counts['ALT'] += alt

                pop_counts[pop] = allele_counts.copy()

            results.append(pop_counts.copy())

        return results

    def slice_vcf(self, vcf_bgzipped_file, chrm, start, stop):

        tbx = pysam.Tabixfile(vcf_bgzipped_file)

        try:
            vcf_slice = tbx.fetch(chrm, start, stop)
        except ValueError:
            return ()

        else:
            chunk = tuple(self.parse_vcf_line(row, self.empty_vcf_line) for row in vcf_slice)
            return ((chrm, start, stop), chunk)
            #return tuple(row for row in vcf_slice)

    def vcf_file_iterator(self, as_dict=True):

        for line in self.__open_vcf__():
            if line.startswith("#") is not True:
                if as_dict == True:
                    yield self.parse_vcf_line(line, self.empty_vcf_line)
                else:
                    yield line
            else:
                continue

    def parse_info_field(self, info_field):

        info_dict = {}
        for item in info_field.split(';'):
            pair = item.split("=")
            if len(pair) == 2:
                info_dict[pair[0]] = pair[1]   # this could be improved on

        return info_dict

    def get_population_sizes(self, vcfline, populations):

        sample_counts = {}

        for pop in populations.keys():
            sample_count = 0

            for sample in populations[pop]:
                if vcfline[sample] is not None:
                    sample_count += 1

            sample_counts[pop] = sample_count

        return sample_counts


    def parse_vcf_line(self, pos, vcf_line_dict):
        """Read in VCF line and convert it to an OrderedDict"""

        pos_parts = pos.strip().split()
        for count, item in enumerate(vcf_line_dict):
            vcf_line_dict[item] = pos_parts[count]

        sample_format = vcf_line_dict["FORMAT"].split(":")

        for count, item in enumerate(vcf_line_dict):

            if count >= 9:
                genotype = vcf_line_dict[item]

                if  "./." in genotype or genotype == ".":      # "'./.'' for dip, '.' for haploid
                    vcf_line_dict[item] = None

                else:
                    genotype = dict(zip(sample_format, genotype.split(":")))

                    # CONVERT STRINGS TO APPOPRIATE TYPES (INTS, FLOATS, ETC.)
                    if genotype.has_key("GQ"):
                        genotype['GQ'] = float(genotype['GQ'])

                    if genotype.has_key("DP"):
                        genotype['DP'] = int(genotype['DP'])

                    if genotype.has_key("AD"):
                        genotype['AD'] = tuple(int(ad) for ad in genotype['AD'].split(","))

                    if genotype.has_key("PL"):
                        genotype['PL'] = tuple(int(ad) for ad in genotype['PL'].split(","))

                    vcf_line_dict[item] = genotype

        vcf_line_dict['POS'] = int(vcf_line_dict['POS'])

        try:
            vcf_line_dict['QUAL'] = float(vcf_line_dict['QUAL'])
        except ValueError:
            pass

        vcf_line_dict['INFO'] = self.parse_info_field(vcf_line_dict['INFO'])

        return vcf_line_dict.copy()


    def filter_vcf_line(self, vcf_line_dict, filters=None):
        parse_vcf_line
        pass


    def lines_2_dicts(self, chunk):

        vcf_line_dict = self.empty_vcf_line.copy()
        return [self.parse_vcf_line(line, vcf_line_dict) for line in chunk]

    def get_chromosome_lengths(self, regions_to_skip=[]):

        try:
            tbx = pysam.Tabixfile(self.input)  # TODO: create try statement to test that file is actually a VCF
        except:
            print 'Input not Tabix Indexed.'
            sys.exit()

        # PARSE LENGTH INFO FROM HEADER

        chrm_lengths = []
        chrm_lengths_dict = {}
        for line in tbx.header:

            if line.startswith("##contig="):

                chrm_name = re.findall(r'ID=.*,', line)
                chrm_name = chrm_name[0].strip('ID=').strip(',')

                chrm_length = re.findall(r'length=.*>', line)
                chrm_length = int(chrm_length[0].strip('length=').strip('>'))

                if chrm_name in regions_to_skip:
                    print 'skipping', chrm_name
                    continue

                chrm_lengths.append((chrm_name, 1, chrm_length))
                chrm_lengths_dict[chrm_name] = chrm_length

        chrm_lengths = tuple(chrm_lengths)
        tbx.close()

        return chrm_lengths_dict

    def get_slice_indicies(self):
        """Get slice information from VCF file that is tabix indexed (bgzipped). """

        # GENERATE SLICES:
        # Does not yield final partial slice. Not a bug!

        #print [c for c in self.chrms_2_sizes.iteritems()]

        if self.region == [None]:

            # ITERATE OVER CHROMOSOMES (USE ORDERED DICT TO KEEP IN VCF ORDER)
            for chrm, length in self.chrms_2_sizes.iteritems():


                cStart = 0
                cStop = 0
                iCount = 0

                while cStop < length:

                    if iCount == 0:
                        cStart = 1
                        cStop = self.window_size
                        iCount += 1

                    yield (chrm, cStart, cStop)
                    cStart += self.step
                    cStop += self.step

        else:

            chrm, start, stop = self.region

            cStart = 0
            cStop = 0
            iCount = 0

            if self.window_size == None:
                self.window_size = stop - start

            if self.step == None:
                self.step = 0

            while cStop < stop:

                if iCount == 0:
                    cStart = start
                    cStop = start + self.window_size - 1
                    iCount += 1

                yield (chrm, cStart, cStop)

                if self.step == 0:
                    cStart += self.window_size
                else:
                    cStart += self.step

                cStop = cStart + self.window_size + self.step - 1

    def snp_chunks_interator(self, num_snv_per_chunk=3, fillvalue=None):

        def grouper(iterable, n):
            """Yields chunks of SNV of n number of SNVs.
               Properly handles chromsome so that
               chunks don't overlap chromosomes."""

            chrom = None
            chunk = []

            for count, i in enumerate(iterable):

                #s = i.split("\t")[:3]
                current_chrom = i["CHROM"]

                if count == 0:
                    chrom = current_chrom

                if current_chrom != chrom:

                    yield chunk
                    chunk = [i]
                    chrom = current_chrom
                    continue

                if current_chrom == chrom:

                    chunk.append(i)

                    if len(chunk) >= n:
                        yield chunk
                        chunk = []

                chrom = current_chrom

        chunk_iterator = grouper(self.vcf_file_iterator(), n=num_snv_per_chunk)

        for c, i in enumerate(chunk_iterator):

            if len(i) == 0: # handles edge case
                continue

            if len(i) != num_snv_per_chunk:
                padding =  num_snv_per_chunk - len(i)
                i.extend([fillvalue]*padding)

            yield i

    def calc_allele_counts(self, vcf_line_dict, sample_ids=None):

        #allele_counts = defaultdict({0:0.0,1:0.0,2:0.0,3:0.0,4:0.0})
        allele_counts = dict((key, {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0}) for key in self.populations.keys())

        for population in self.populations.keys():

            for sample_id in self.populations[population]:

                if vcf_line_dict[sample_id] != None:

                    genotype = vcf_line_dict[sample_id]
                    genotype = genotype["GT"].split("/")   # TODO add phased logic

                    if genotype == [".", "."]:
                        continue

                    genotype = [int(item) for item in genotype]

                    for allele in genotype:
                        allele_counts[population][allele] += 1.0

        return allele_counts







