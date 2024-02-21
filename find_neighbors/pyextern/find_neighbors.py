import ctypes
import os


class FindNeighbors:

    def __init__(self,
                 so_path: str = None):
        # TODO prepare a default path for so file
        if so_path is None:
            so_path = "find_neighbors_python_external.so"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.handle = ctypes.CDLL(os.path.join(dir_path, so_path))
        self.handle.find_neighbors.argtypes = [ctypes.c_int,
                                               ctypes.c_int,
                                               ctypes.c_float,
                                               ctypes.c_char_p,
                                               ctypes.c_char_p,
                                               ]

    def run(self,
            batch_index: int,
            n_batch: int,
            z_max: float,
            input_file: str,
            output: str):
        """The string input to char* should be in bytes, encode the byte result"""
        return self.handle.find_neighbors(batch_index,
                                          n_batch,
                                          z_max,
                                          input_file.encode(),
                                          output.encode())


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="find_neighbors.py",
                                     description="find neighbors for normalized mosdepth result",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-C', "--index_count", type=int, required=True,
        help="the index of batch (0-based), for parallel setting use")
    parser.add_argument(
        '-N', "--n_batch", type=int, required=True,
        help="total number of batches included")
    parser.add_argument(
        '-Z', "--z_max", type=float, required=True,
        help="the maximum range for z-score, any value exceed will be cropped to z_max")
    parser.add_argument(
        "-I", "--input_file", type=str, required=True,
        help="the path of the input file, should be a *_ID_scale_zdepths.txt.gz file")
    parser.add_argument(
        "-O", "--output_prefix", type=str, required=True,
        help="the output path prefix, the output will be <output_prefix>_zMax_<z_max>.txt.gz")
    Args = parser.parse_args()
    fn = FindNeighbors()
    fn.run(batch_index=Args.index_count,
           n_batch=Args.n_batch,
           z_max=Args.z_max,
           input_file=Args.input_file,
           output=Args.output_prefix
           )
