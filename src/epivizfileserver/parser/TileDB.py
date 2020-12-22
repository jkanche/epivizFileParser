import tiledb
import numpy as np
import pandas as pd
import ujson
from .TileDBTbxFile import TileDBTbxFile

class TileDB(object):
    """
    TileDB Class to parse only local tiledb files 

    Args:
        path (str): local full path to a dataset tiledb_folder. This folder
            should contain data.tiledb, rows and cols files. See below for more detail.
        columns ([str]) : column names for various columns in file

    Detail:
        The tiledb_folder should contain:
            'data.tiledb' directory - corresponds to the uri of a tiledb array. The tiledb array
            must have a 'vals' attribute from which values are read. The array should have as many
            rows as the number of lines in the 'rows' file, and as many columns as the number of
            lines in the 'cols' file.

            'rows' file - this is a tab-separated value file describing the rows of the tiledb array
            it must have as many lines as rows in the tiledb file. There should be no index column in
            this file (i.e., it is read with pandas.read_csv(..., sep='\t', index_col=False)). It must
            have columns 'chr', 'start' and 'end'.
            We index the rows file using Tabix so we are not loading the entire file into memory.
            This file contains columns as annotated in .json file

            'cols' file - this is a tab-separated value file describing the columns of the tiledb array.
            It must have as many files as columns in the tiledb file. Column names for the tiledb array
            will be obtained from the first column in this file (i.e., iti is read with 
            pandas.read_csv(..., sep='\t', index_col=0)). 
    """
    def __init__(self, path):
        self.path = path
        self.count = tiledb.open(path + "/data.tiledb", 'r')

        # get columns of the rows file
        row_rec = ujson.load(open(path + "/rows.tsv.bgz.json"))
        metadata = [m["name"] for m in row_rec["covariates"]]
        fmeta = []

        # renaming columns
        for m in metadata:
            if m.lower() == "id":
                m = "gene"

            if m.lower() == "seqnames":
                m = "chr"
            fmeta.append(m)

        # metadata = [m for m in metadata if m not in ['seqnames', 'start', 'end', 'chr']]
        self.rows = TileDBTbxFile(path + "/rows.tsv.bgz", columns=fmeta)
        self.cols = pd.read_csv(path + "/cols.tsv", sep="\t", index_col=0)
        self.columns = self.cols["epiviz_ids"].values # self.cols.index.values

    def getRange(self, chr, start = None, end = None, bins=2000, zoomlvl=-1, metric="AVG", respType = "DataFrame", treedisk=None):
        """Get data for a given genomic location

        Args:
            chr (str): chromosome 
            start (int): genomic start
            end (int): genomic end
            respType (str): result format type, default is "DataFrame

        Returns:
            result
                a DataFrame with matched regions from the input genomic location if respType is DataFrame else result is an array
            error 
                if there was any error during the process
        """
        result = pd.DataFrame(columns=self.columns)

        try:
            # result_rows = self.rows[(self.rows["chr"] == chr) & (self.rows["start"] <= end) & (self.rows["end"] >= start)]
            result_rows, err = self.rows.getRange('"' + chr + '"', start, end)
            result_rows = result_rows.applymap(lambda x: x.replace('"', ''))
            
            indices = result_rows["X__rowindex"].values.astype(int)
            matrix = self.count[min(indices):max(indices)+1,]['vals']
            
            result_matrix = pd.DataFrame(matrix, index=range(min(indices), max(indices)+1), columns=self.columns)
            result_merge = pd.concat([result_rows, result_matrix], axis=1, join="inner")
            return result_merge, None
        except Exception as e:
            print(str(e))
            return result, str(e)
