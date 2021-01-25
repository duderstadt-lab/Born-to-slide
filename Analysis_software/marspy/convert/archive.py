import pandas as pd
from jnius import autoclass

from awesome_data import DataSet
from marspy.convert.molecule import *


class Archive:

    def __init__(self, filepath):
        self.filepath = filepath
        self.name = self.filepath.split('/')[-1]
        self.File = autoclass('java.io.File')
        self.yamaFile = self.File(self.filepath)

    def get_molecule_by_uid(self, uid):
        raise NotImplementedError

    def get_molecules_by_tags(self, tags):
        raise NotImplementedError

    def validate_params(self):
        pass


class SingleMoleculeArchive(Archive):
    instances = []

    def __init__(self, filepath, accept_tag, label=dict()):
        Archive.__init__(self, filepath)
        self.instances.append(self)
        self.Archive = autoclass('de.mpg.biochem.mars.molecule.SingleMoleculeArchive')
        self.archive_link = self.Archive(self.yamaFile)
        self.metadata_uids = list(self.archive_link.getMetadataUIDs())
        self.label = label

        # nucleotide
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('nucleotide')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nucleotide')) == 0:
            # default n/a
            self.nucleotide = 'n/a'
            print(f'nucleotide not found. Setting default to {self.nucleotide}')
        # parameter properly set
        else:
            self.nucleotide = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nucleotide')

        # highsalt_wash
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getParameter('highsalt_wash')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if Parameter is not set, getParameter returns np.nan
        if np.isnan(self.archive_link.getMetadata(self.metadata_uids[0]).getParameter('highsalt_wash')):
            # default False
            self.highsalt_wash = False
            print(f'highsalt_wash not found. Setting default to {self.highsalt_wash}')
        else:
            self.highsalt_wash = \
                self.archive_link.getMetadata(self.metadata_uids[0]).getParameter('highsalt_wash') == 1

        # cdc6
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('cdc6')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('cdc6')) == 0:
            # default n/a
            self.cdc6 = 'n/a'
            print(f'cdc6 not found. Setting default to {self.cdc6}')
        # parameter properly set
        else:
            self.cdc6 = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('cdc6')

        self.protein = list(self.label.keys())[0]

        # instantiate a new SingleMolecule for each uid and store instances as list
        self.molecules = [SingleMolecule(uid, self.protein, archive=self.archive_link) for uid in
                          self.archive_link.getMoleculeUIDs() if self.archive_link.get(uid).hasTag(accept_tag)]
        self.tags = set()
        for molecule in self.molecules:
            self.tags.update(molecule.tags)

    def get_molecule_by_uid(self, uid):
        """
        Returns molecule object with provided UID.
        """
        return list(filter(lambda molecule: molecule.uid == uid, self.molecules))[0]

    def get_molecules_by_tags(self, tags):
        """
        Provide tags as list.
        Returns filter of all molecules which have all the specified tags
        """
        return filter(lambda molecule: set(tags).issubset(set(molecule.tags)), self.molecules)

    def __len__(self):
        return len(self.molecules)


class DnaMoleculeArchive(Archive):
    instances = []

    def __init__(self, filepath, accept_tag, labels=dict()):
        Archive.__init__(self, filepath)
        self.instances.append(self)
        self.Archive = autoclass('de.mpg.biochem.mars.molecule.DnaMoleculeArchive')
        self.archive_link = self.Archive(self.yamaFile)
        self.metadata_uids = list(self.archive_link.getMetadataUIDs())
        self.dna_molecule_count = 0
        for metadata in self.metadata_uids:
            self.dna_molecule_count += dict(sc.to_python(self.archive_link.getMetadata(metadata).getParameters()))[
                'DnaMoleculeCount']
        # subtract # of reject_dna tags
        self.dna_molecule_count -= len(list(filter(lambda uid:
                                                   self.archive_link.get(uid).hasTag('reject_dna'),
                                                   self.archive_link.moleculeUIDs)))
        self.labels = labels
        # nucleotide
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('nucleotide')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nucleotide')) == 0:
            # default n/a
            self.nucleotide = 'n/a'
            print(f'nucleotide not found. Setting default to {self.nucleotide}')
        # parameter properly set
        else:
            self.nucleotide = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nucleotide')

        # highsalt_wash
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getParameter('highsalt_wash')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if Parameter is not set, getParameter returns np.nan
        if np.isnan(self.archive_link.getMetadata(self.metadata_uids[0]).getParameter('highsalt_wash')):
            # default False
            self.highsalt_wash = False
            print(f'highsalt_wash not found. Setting default to {self.highsalt_wash}')
        else:
            self.highsalt_wash = \
                self.archive_link.getMetadata(self.metadata_uids[0]).getParameter('highsalt_wash') == 1

        # dna_count_valid: data was fully analyzed - ALL DNA molecules fitted
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getParameter('dna_count_valid')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if Parameter is not set, getParameter returns np.nan
        if np.isnan(self.archive_link.getMetadata(self.metadata_uids[0]).getParameter('dna_count_valid')):
            # default True
            self.dna_count_valid = True
            print(f'dna_count_valid not found. Setting default to {self.dna_count_valid}')
        else:
            self.dna_count_valid = \
                self.archive_link.getMetadata(self.metadata_uids[0]).getParameter('dna_count_valid') == 1

        # t7_terminator
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getParameter('t7_terminator')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if Parameter is not set, getParameter returns np.nan
        if np.isnan(self.archive_link.getMetadata(self.metadata_uids[0]).getParameter('t7_terminator')):
            # default False
            self.t7_terminator = False
            print(f't7_terminator not found. Setting default to {self.t7_terminator}')
        else:
            self.t7_terminator = \
                self.archive_link.getMetadata(self.metadata_uids[0]).getParameter('t7_terminator') == 1

        # chromatin
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('chromatin')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('chromatin')) == 0:
            # default n/a
            self.chromatin = 'n/a'
            print(f'chromatin not found. Setting default to {self.chromatin}')
        # parameter properly set
        else:
            self.chromatin = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('chromatin')

        # cdc6
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('cdc6')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('cdc6')) == 0:
            # default n/a
            self.cdc6 = 'n/a'
            print(f'cdc6 not found. Setting default to {self.cdc6}')
        # parameter properly set
        else:
            self.cdc6 = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('cdc6')

        self.proteins = set()
        # will find all columns in metadata DataTable with 'ChannelIndex [Protein]'
        for match in re.findall('ChannelIndex \w+', '$'.join(dict(sc.to_python(
                self.archive_link.getMetadata(self.metadata_uids[0]).getImage(0).getPlane(0, 0,
                                                                                          0).getStringFields())).keys())):
            self.proteins.add(match.split()[-1])

        # instantiate a new DnaMolecule for each uid and store instances as list
        self.molecules = [DnaMolecule(uid, self.proteins, archive=self.archive_link) for uid in
                          self.archive_link.getMoleculeUIDs()
                          if self.archive_link.get(uid).hasTag(accept_tag)]

        # define archive tags union of all molecule tags
        # define archive prefixes as union of all molecule prefixes (will be used for top level columns in big df later)
        self.tags = set()
        self.prefixes = set()
        for molecule in self.molecules:
            self.tags.update(molecule.tags)
        self.prefixes.update(molecule.prefixes)

    def validate_params(self):
        """
        Integrity check of passed Archive.
        """
        # compare number protein in params vs actual one (retrieved from metadata)
        for molecule in self.molecules:
            # take global protein to confirm dict was pasted correctly
            for protein in self.proteins:
                if not (molecule.proteins[protein] == molecule.params['Number_' + protein]):
                    err_message = f"Conflict in molecule {molecule.uid}!\n\
                    Number of {protein} retrieved from metadata: {molecule.proteins[protein]}\n\
                    Number of {protein} based on Parameter: {molecule.params['Number_' + protein]}"
                    raise MarsPyException(err_message)

        return 'passed'

    def add_segments_tables(self):
        """
        Attach all segment tables to molecule records (stored as dict)
        """
        # Do we have collisions in the archive?
        coll_exp = False
        for tag in self.tags:
            if re.match('coll', tag):
                coll_exp = True

        for molecule in self.molecules:
            molecule.seg_dfs = list()

            # all segmentTableNames
            for x, y, region in (sc.to_python(self.archive_link.get(molecule.uid).getSegmentsTableNames())):
                # internal control that all seg_dfs are valid
                _assigned = False
                # all proteins on molecule
                for prefix in molecule.prefixes:
                    if re.match(prefix, x) and re.match(prefix, y):
                        molecule.seg_dfs.append(SegmentsTable(molecule=molecule, prefix=prefix,
                                                              col_x=x, col_y=y, region=region, coll_exp=coll_exp))
                        _assigned = True
                        break
                if not _assigned:
                    err_message = f"Conflict in molecule {molecule.uid}!\nSegmentTable {x} {y} {region} not assigned!"
                    raise MarsPyException(err_message)

    def filter_segments(self, b_min=0, sigma_b_max=0):
        """
        Filter all segments for all molecules in archive based on SegmentsTable type.
        Also see filter_segments() in SegmentsTable object:

            Mode 1: SegmentTable type: 'bleaching' -
            Reject all steps with increase in fluorescence intensity (initial quenching).
            If increase in fluorescence is detected, remove segments with same intensity (double counts).

            Mode 2: SegmentsTable type: 'rate' -
            Reject all segments with B value (velocity) < b_min and sigma_B > sigma_b_max (poor fits)

        """
        for molecule in self.molecules:
            for prefix in molecule.prefixes:
                # check if one protein molecule has only one seg_df with type 'bleaching'
                if len(list(filter(lambda df:
                                   prefix == df.prefix and df.type == 'bleaching', molecule.seg_dfs))) > 1:
                    err_message = f"Conflict in molecule {molecule.uid}!\nMore than one SegmentTable for {prefix}."
                    raise MarsPyException(err_message)

            # apply filter to all seg_dfs
            for seg_df in molecule.seg_dfs:
                seg_df.filter_segments(b_min=b_min, sigma_b_max=sigma_b_max)

            # in case seg_df is empty after filtering, delete object
            remove_seg_dfs = set()
            for seg_df in molecule.seg_dfs:
                if len(seg_df.df) == 0:
                    remove_seg_dfs.add(seg_df)
            for seg_df in remove_seg_dfs:
                molecule.seg_dfs.remove(seg_df)

    def calc_bleaching_steps(self):
        """
        Calculate bleaching steps for all proteins of all molecules in archive.
        Also see calc_bleaching_steps() in SegmentsTable object:

            Calculate number of bleaching steps based off of segment table rows.
            Returns: number of bleaching steps (integer value)

        No return value here; bleaching steps are stored as attribute dict with prefixes as key.
        """
        for molecule in self.molecules:
            molecule.bleaching_steps = dict()
            for prefix in molecule.prefixes:
                # check if one protein molecule has only one seg_df with type 'bleaching'
                if len(list(filter(lambda seg_df:
                                   prefix == seg_df.prefix and seg_df.type == 'bleaching', molecule.seg_dfs))) > 1:
                    err_message = f"Conflict in molecule {molecule.uid}!\nMore than one SegmentTable for {prefix}."
                    raise MarsPyException(err_message)

                # only molecules with proper bleaching (not rejected)
                if 'reject_bleach_' + prefix in molecule.tags:
                    continue

                molecule.bleaching_steps[prefix] = \
                    list(filter(lambda seg_df: prefix == seg_df.prefix and seg_df.type == 'bleaching',
                                molecule.seg_dfs))[0].calc_bleaching_steps()

    def detect_pauses(self, thresh=3, global_thresh=False, col='B'):
        """
        Detect pauses in translocation for all SegmentTables of all molecules in archive.
        Also see detect_pauses() in SegmentsTable object:

            Detection pauses in SegmentTable (only for type = 'rate', others are skipped)
            global_thresh: Set to True if a fixed threshold for all molecules should be used
            thresh: threshold to detect pauses.
            If global_thresh is False, a molecule-specific threshold is calculated with thresh^-1 * np.mean(col)
            col: column evaluated for pauses
        """
        for molecule in self.molecules:
            for seg_df in molecule.seg_dfs:
                seg_df.detect_pauses(thresh=thresh, global_thresh=global_thresh, col=col)

    def get_molecule_by_uid(self, uid):
        """
        Returns molecule object with provided UID.
        """
        return list(filter(lambda molecule: molecule.uid == uid, self.molecules))[0]

    def get_molecules_by_tags(self, tags):
        """
        Provide tags as list.
        Returns filter of all molecules which have all the specified tags
        """
        return filter(lambda molecule: set(tags).issubset(set(molecule.tags)), self.molecules)

    def __len__(self):
        return len(self.molecules)


def instantiate_archive(name, datasets):
    """
    Instantiates passed archive from underlying dataset
    """
    # check if we have the right data type
    for data in datasets:
        if not isinstance(data, DataSet):
            raise MarsPyException('Dataset contains non-compatible data type.')
    data = list(filter(lambda dataset: dataset.name == name, datasets))[0]
    if data.archive_type == 'DnaMoleculeArchive':
        DnaMoleculeArchive(filepath=data.filepath + data.name, accept_tag=data.accept_tag, labels=data.labels)
    elif data.archive_type == 'SingleMoleculeArchive':
        SingleMoleculeArchive(filepath=data.filepath + data.name, accept_tag=data.accept_tag, label=data.labels)
    else:
        raise MarsPyException(f'Failed to instantiate Archive {data.name}.')


def describe_archives(archives):
    """
    Describes passed archives by returning a pandsa DataFrame. Pass archives as iterable object
    """
    df = pd.DataFrame(columns=['# of datasets', '# of molecules', 'labeled proteins', 'nucleotide', 'HS challenge?',
                               'chromatin', 'terminator?', 'archive validation'])
    for archive in archives:
        _temp_df = pd.DataFrame(index=[archive.name.split('.')[0]],
                                data=[[len(archive.metadata_uids), len(archive),
                                       '; '.join([label + '-' + protein for protein, label in
                                                  DnaMoleculeArchive.instances[0].labels.items()]),
                                       archive.nucleotide, archive.highsalt_wash, archive.chromatin,
                                       archive.t7_terminator, archive.validate_params()]],
                                columns=['# of datasets', '# of molecules', 'labeled proteins', 'nucleotide',
                                         'HS challenge?', 'chromatin', 'terminator?', 'archive validation'])
        df = pd.concat([df, _temp_df])
    df = df.infer_objects()
    return df
