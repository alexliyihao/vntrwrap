import ctypes
import os


class NormalizeMosdepth:

    def __init__(self,
                 so_path: str = None):
        # TODO prepare a default path for so file
        if so_path is None:
            so_path = "normalize_mosdepth_python_external.so"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.handle = ctypes.CDLL(os.path.join(dir_path, so_path))
        self.handle.normalize_mosdepth.argtypes = [ctypes.c_char_p,
                                                   ctypes.c_char_p,
                                                   ctypes.c_char_p,
                                                   ctypes.c_char_p,
                                                   ctypes.c_int]

    def run(self,
            mosdepth_prefix: str,
            ucsc_bed_file: str,
            example_input: str,
            output_path: str,
            n_sample: int):
        """The string input to char* should be in bytes, encode the byte result"""
        return self.handle.normalize_mosdepth(mosdepth_prefix.encode(),
                                              ucsc_bed_file.encode(),
                                              example_input.encode(),
                                              output_path.encode(),
                                              n_sample)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="normalize_mosdepth.py",
                                     description="normalize the output of mosdepth analysis",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-I', "--mosdepth_prefix", type=str, required=True,
        help="prefix of mosdepth input (no more than 170 characters), used in <prefix>_batch_<batchnumber>.txt.gz")
    parser.add_argument(
        '-B', "--bed_file", type=str, required=True,
        help="bed file path e.g. /path/to/repeat_mask_list.hg38.ucsc_bed")
    parser.add_argument(
        '-O', "--output_path", type=str, required=True,
        help="output path e.g. /path/to/ID_scale_zdepths.txt.gz")
    parser.add_argument(
        "-E", "--example", type=str, required=True,
        help="example mosdepth output as input e.g. /path/to/name_regions.bed.gz")
    parser.add_argument(
        "-N", "--n_sample", type=int, required=True,
        help="the sample size, in integer")
    Args = parser.parse_args()
    nm = NormalizeMosdepth()
    nm.run(mosdepth_prefix=Args.mosdepth_prefix,
           ucsc_bed_file=Args.bed_file,
           example_input=Args.example,
           output_path=Args.output_path,
           n_sample=Args.n_sample
           )
