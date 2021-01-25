import re
from collections import defaultdict

import numpy as np
import scyjava as sc
import seaborn as sns
from scyjava.convert._pandas import table_to_pandas


class Molecule:

    def __init__(self, uid, archive):
        self.uid = uid
        self.archive = archive
        self.meta_uid = self.archive.get(self.uid).getMetadataUID()
        self.params = dict(sc.to_python(self.archive.get(self.uid).getParameters()))
        self.tags = list(sc.to_python(self.archive.get(self.uid).getTags()))
        self.regions = list(sc.to_python(self.archive.get(self.uid).getRegionNames()))
        # don't keep entire df in memory for future versions
        self.df = table_to_pandas(self.archive.get(self.uid).getTable())
        self.seg_dfs = None

    def __str__(self):
        return f'Greetings from Molecule {self.uid}.'

    def __len__(self):
        return len(self.df)

    def __del__(self):
        pass
        print(f'Molecule {self.uid} deleted.')


class SingleMolecule(Molecule):

    def __init__(self, uid, protein, archive):
        Molecule.__init__(self, uid, archive)
        self.protein = protein


class DnaMolecule(Molecule):

    def __init__(self, uid, proteins, archive):
        Molecule.__init__(self, uid, archive)

        # DnaMolecule specific attributes
        self.proteins = {protein: 0 for protein in proteins}

        # grab keys from dict
        for protein in self.proteins:
            # protein specific prefixes with nomenclature protein_prefixes:
            (exec(
                f"self.{protein + '_prefixes'} = set(re.findall(f'{protein}_\d+_',' '.join(self.df.columns)))"))

            # Store number of molecules based off of actual dataTable headers
            self.proteins[protein] = len(set(re.findall(f'{protein}_\d+_', ' '.join(self.df.columns))))

        # generate prefixes based union of protein_prefixes
        self.prefixes = set()
        for protein in self.proteins.keys():
            for prefix in getattr(self, f'{protein}_prefixes'):
                self.prefixes.add(prefix)

        # region objects for DnaMolecules
        self.regions = list()
        # all region names
        for region_name in sc.to_python(self.archive.get(self.uid).getRegionNames()):
            _region = self.archive.get(self.uid).getRegion(region_name)
            match_prefix = None
            # separate prefix from column name
            for prefix in self.prefixes:
                if re.match(prefix, _region.getColumn()):
                    # correct prefix found
                    match_prefix = prefix
                    break

            # append Region object to molecule regions
            self.regions.append(Region(uid=self.uid,
                                       name=region_name,
                                       start=_region.getStart(),
                                       end=_region.getEnd(),
                                       prefix=match_prefix,
                                       column=_region.getColumn().split(match_prefix)[-1]))

    def calc_length_dna(self):
        """
        Calculates the Molecule's DNA length in px.
        """
        return np.sqrt((self.params['Dna_Bottom_x2'] - self.params['Dna_Top_x1']) ** 2 +
                       (self.params['Dna_Bottom_y2'] - self.params['Dna_Top_y1']) ** 2)

    def plot(self):
        for prefix in self.prefixes:
            try:
                sns.lineplot(x=prefix + 'Time (s)', y=prefix + 'Position_on_DNA', data=self.df)
            except ValueError:
                continue


class SegmentsTable:
    """
    SegmentsTable object holding the actual df, with additional information in attributes.
    Also contains specific methods for filtering, bleaching steps, pause detection.
    """

    def __init__(self, molecule, prefix, col_x, col_y, region, coll_exp=False):
        # uid for debugging
        self.uid = molecule.uid
        # which prefix does it belong to
        self.prefix = prefix
        # remove prefix for x and y columns => general naming
        self.col_x = col_x.split(self.prefix)[-1]
        self.col_y = col_y.split(self.prefix)[-1]
        self.region = region
        # actual SegmentsTable()
        self.df = sc.to_python(molecule.archive.get(molecule.uid).getSegmentsTable(col_x, col_y, region))
        # multicolor experiments with transcription
        self.coll_exp = coll_exp
        # type of SegmentTable (default False)
        self.type = None
        # keep track if seg_dfs were already filtered before data is interpreted
        self.filtered = False
        # assign type of SegmentTable (one needs to be True)
        if self.col_y == 'Intensity':
            self.type = 'bleaching'
        elif self.col_y == 'Position_on_DNA':
            self.type = 'rate'
        else:
            err_message = f"Conflict in molecule {self.uid}!\nSegmentTable type could not be assigned!"
            raise MarsPyException(err_message)

        # Information about proteins pushed is in the region name with format: *number**protein1*-
        # default: empty dictionary
        self.pushed_proteins = defaultdict(int)
        if self.coll_exp:
            # region nomenclature: region_*PREFIX*_*protein1*:n-*protein2*:n; e.g.region_T7_1_NUC:2-MCM:1
            for pushed_protein in self.region.split(self.prefix)[-1].split('-'):
                try:
                    self.pushed_proteins[pushed_protein.split(':')[0]] = int(pushed_protein.split(':')[1])
                # no collision / pushing for this SegmentsTable instance
                except IndexError:
                    break

    def filter_segments(self, b_min=0, sigma_b_max=0):
        """
        Mode 1: SegmentTable type: 'bleaching' -
        Reject all steps with increase in fluorescence intensity (initial quenching).
        If increase in fluorescence is detected, remove segments with same intensity (double counts).

        Mode 2: SegmentsTable type: 'rate' -
        Reject all segments with B value (velocity) < b_min and sigma_B > sigma_b_max (poor fits)
        """

        # Mode 1 - type 'bleaching'
        if self.type == 'bleaching':
            step_increase = False

            for i in range(1, len(self.df)):
                if self.df.loc[i, 'y1'] > self.df.loc[i - 1, 'y1']:
                    # set to True to trigger subsequent filtering
                    step_increase = True
                    self.df.drop(i - 1, axis=0, inplace=True)
            # reset index to start with 0 again
            # reset index seg_df
            self.df.reset_index(drop=True, inplace=True)

            if step_increase:
                # need to remove rows at the end
                remove_rows = set()
                # loop through rows and find segments with same intensity (within 4E3)
                for i in range(len(self.df)):
                    # all subsequent rows
                    for j in range(i + 1, len(self.df)):
                        # compare the two rows
                        if abs(self.df.loc[i, 'y1'] - self.df.loc[j, 'y1']) < 4000:
                            remove_rows.add(j)

                # remove all rows & update indices
                self.df.drop(list(remove_rows), axis=0, inplace=True)
                self.df.reset_index(drop=True, inplace=True)

            # successfully filtered
            self.filtered = True

        # Mode 2 - type 'rate'
        elif self.type == 'rate':

            # need to remove rows at the end
            remove_rows = set()
            # loop through rows and find segments which match exclusion criteria
            for i in range(len(self.df)):
                if (self.df.loc[i, 'B'] <= b_min) or (self.df.loc[i, 'sigma_B'] >= sigma_b_max):
                    remove_rows.add(i)

            # remove all rows & update indices
            self.df.drop(list(remove_rows), axis=0, inplace=True)
            self.df.reset_index(drop=True, inplace=True)

            # successfully filtered
            self.filtered = True

    def calc_bleaching_steps(self):
        """
        Calculate number of bleaching steps based off of segment table rows.
        Returns: number of bleaching steps (integer value)
        """
        # check if bleaching analysis is appropriate
        if self.type != 'bleaching':
            err_message = f"Conflict in molecule {self.uid}!\n\
            A SegmentTable type is {self.type} but calculation of bleaching steps was requested!"
            raise MarsPyException(err_message)

        # Filter already applied?
        elif not self.filtered:
            print(f"Warning! A SegmentsTable in molecule {self.uid} was analyzed without filters applied!")

        return len(self.df) - 1

    def detect_pauses(self, thresh=3, global_thresh=False, col='B'):
        """
        Detection pauses in SegmentTable (only for type = 'rate')
        global_thresh: Set to True if a fixed threshold for all molecules should be used
        thresh: threshold to detect pauses.
        If global_thresh is False, a molecule-specific threshold is calculated with thresh^-1 * np.mean(col)
        col: column evaluated for pauses
        """
        # only SegmentsTables with type 'rate'
        if self.type == 'rate':
            self.df['pause_' + col] = False
            cutoff = thresh
            # iterate once for each entry in seg_df
            for i in range(len(self.df)):
                if not global_thresh:
                    # redefine cutoff
                    cutoff = self.df[~self.df['pause_' + col]][col].mean() / thresh

                for row in self.df.index:
                    self.df.loc[row, 'pause_' + col] = self.df.loc[row, col] < cutoff
            # if two subsequent segments are pauses merge them
            remove_rows = set()
            for i in range(1, len(self.df)):
                # both pauses
                # additional requirements: time values match and y values within 1 kb
                if (self.df.loc[i - 1, 'pause_B'] and self.df.loc[i, 'pause_B'] and
                        self.df.loc[i - 1, 'x2'] == self.df.loc[i, 'x1'] and
                        abs(self.df.loc[i - 1, 'y2'] - self.df.loc[i, 'y1'] < 1000)):
                    # recalculate values (x2 and y2 values in row i stay the same)
                    x1 = self.df.loc[i - 1, 'x1']
                    y1 = self.df.loc[i - 1, 'y1']
                    a = np.average(self.df[i - 1:i + 1]['A'],
                                   weights=self.df[i - 1:i + 1]['x2'] - self.df[i - 1:i + 1]['x1'])
                    sigma_a = np.average(self.df[i - 1:i + 1]['sigma_A'],
                                         weights=self.df[i - 1:i + 1]['x2'] - self.df[i - 1:i + 1]['x1'])
                    b = np.average(self.df[i - 1:i + 1]['B'],
                                   weights=self.df[i - 1:i + 1]['x2'] - self.df[i - 1:i + 1]['x1'])
                    sigma_b = np.average(self.df[i - 1:i + 1]['sigma_B'],
                                         weights=self.df[i - 1:i + 1]['x2'] - self.df[i - 1:i + 1]['x1'])

                    # add row i-1 for removal and reassign calculated values for row i
                    remove_rows.add(i - 1)
                    self.df.loc[i, 'x1'] = x1
                    self.df.loc[i, 'y1'] = y1
                    self.df.loc[i, 'A'] = a
                    self.df.loc[i, 'sigma_A'] = sigma_a
                    self.df.loc[i, 'B'] = b
                    self.df.loc[i, 'sigma_B'] = sigma_b

            # remove all rows & update indices
            self.df.drop(list(remove_rows), axis=0, inplace=True)
            self.df.reset_index(drop=True, inplace=True)

    def calc_rate(self, burst=True):
        """
        Calculate protein's (burst) velocity based on time-weighted segments.
        Input: burst (default True): if True, ignore pause segments for calculation (requires previous pause detection)
        Returns: rate (weighted average of rate segments), time (total time of used segments)
        """
        # check if rate calculation is appropriate
        if self.type != 'rate':
            err_message = f"Conflict in molecule {self.uid}!\n\
                    A SegmentTable type is {self.type} but rate calculation was requested!"
            raise MarsPyException(err_message)

        # burst mode
        if burst:
            # check if pause detection was run before
            if 'pause_B' not in self.df.columns:
                err_message = f"Conflict in molecule {self.uid}!\n\
                Burst rate calculation requires previous pause detection (see detect_pauses())."
                raise MarsPyException(err_message)

            # Filter already applied?
            elif not self.filtered:
                print(f"Warning! A SegmentsTable in molecule {self.uid} was analyzed without filters applied!")

            rate = np.average(self.df[~self.df['pause_B']]['B'],
                              weights=self.df[~self.df['pause_B']]['x2'] - self.df[~self.df['pause_B']]['x1'])
            time = (self.df[~self.df['pause_B']]['x2'] - self.df[~self.df['pause_B']]['x1']).sum()
            return rate, time

        # non-burst mode
        else:
            if not self.filtered:
                print(f"Warning! A SegmentsTable in molecule {self.uid} was analyzed without filters applied!")

            rate = np.average(self.df['B'], weights=self.df['x2'] - self.df['x1'])
            time = (self.df['x2'] - self.df['x1']).sum()
            return rate, time


class Region:
    """
    Region object holding following attributes: name, start, end and molecule column.
    """

    def __init__(self, uid, name, start, end, prefix, column):
        # uid for debugging
        self.uid = uid
        self.name = name
        self.start = start
        self.end = end
        self.prefix = prefix
        # to which column was the region set (without the prefix)
        self.column = column


# a custom error if something goes wrong
class MarsPyException(Exception):
    def __init___(self):
        Exception.__init__(self)


class MarsPyWarning(UserWarning):
    def __init__(self):
        UserWarning.__init__(self)
