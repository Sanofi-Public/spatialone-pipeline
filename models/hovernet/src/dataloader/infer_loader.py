"""
infer_loader.py
"""
import math
import sys

import cv2
import matplotlib.pyplot as plt
import numpy as np
import psutil
import torch
import torch.utils.data as data


####
class SerializeFileList(data.IterableDataset):
    """Read a single file as multiple patches of same shape, perform the padding beforehand."""

    def __init__(self, img_list, patch_info_list, patch_size, preproc=None):
        """init

        Args:
            img_list (_type_): _description_
            patch_info_list (_type_): _description_
            patch_size (_type_): _description_
            preproc (_type_, optional): _description_. Defaults to None.
        """
        super().__init__()
        self.patch_size = patch_size

        self.img_list = img_list
        self.patch_info_list = patch_info_list

        self.worker_start_img_idx = 0
        # * for internal worker state
        self.curr_img_idx = 0
        self.stop_img_idx = 0
        self.curr_patch_idx = 0
        self.stop_patch_idx = 0
        self.preproc = preproc
        return

    def __iter__(self):
        """iter

        Returns:
            _type_: _description_
        """
        worker_info = torch.utils.data.get_worker_info()
        if worker_info is None:  # single-process data loading, return the full iterator
            self.stop_img_idx = len(self.img_list)
            self.stop_patch_idx = len(self.patch_info_list)
            return self
        else:  # in a worker process so split workload, return a reduced copy of self
            per_worker = len(self.patch_info_list) / float(worker_info.num_workers)
            per_worker = int(math.ceil(per_worker))

            global_curr_patch_idx = worker_info.id * per_worker
            global_stop_patch_idx = global_curr_patch_idx + per_worker
            self.patch_info_list = self.patch_info_list[
                global_curr_patch_idx:global_stop_patch_idx
            ]
            self.curr_patch_idx = 0
            self.stop_patch_idx = len(self.patch_info_list)
            # * check img indexer, implicit protocol in infer.py
            global_curr_img_idx = self.patch_info_list[0][-1]
            global_stop_img_idx = self.patch_info_list[-1][-1] + 1
            self.worker_start_img_idx = global_curr_img_idx
            self.img_list = self.img_list[global_curr_img_idx:global_stop_img_idx]
            self.curr_img_idx = 0
            self.stop_img_idx = len(self.img_list)
            return self  # does it mutate source copy?

    def __next__(self):
        """next

        Raises:
            StopIteration: _description_

        Returns:
            _type_: _description_
        """
        if self.curr_patch_idx >= self.stop_patch_idx:
            raise StopIteration  # when there is nothing more to yield
        patch_info = self.patch_info_list[self.curr_patch_idx]
        img_ptr = self.img_list[patch_info[-1] - self.worker_start_img_idx]
        patch_data = img_ptr[
            patch_info[0] : patch_info[0] + self.patch_size,
            patch_info[1] : patch_info[1] + self.patch_size,
        ]
        self.curr_patch_idx += 1
        if self.preproc is not None:
            patch_data = self.preproc(patch_data)
        return patch_data, patch_info


####
class SerializeArray(data.Dataset):
    """SerializeArray

    Args:
        data (_type_): _description_
    """

    def __init__(self, mmap_array_path, patch_info_list, patch_size, preproc=None):
        """init

        Args:
            mmap_array_path (_type_): _description_
            patch_info_list (_type_): _description_
            patch_size (_type_): _description_
            preproc (_type_, optional): _description_. Defaults to None.
        """
        super().__init__()
        self.patch_size = patch_size

        # use mmap as intermediate sharing, else variable will be duplicated
        # accross torch worker => OOM error, open in read only mode
        self.image = np.load(mmap_array_path, mmap_mode="r")

        self.patch_info_list = patch_info_list
        self.preproc = preproc
        return

    def __len__(self):
        """len

        Returns:
            _type_: _description_
        """
        return len(self.patch_info_list)

    def __getitem__(self, idx):
        """getitem

        Args:
            idx (_type_): _description_

        Returns:
            _type_: _description_
        """
        patch_info = self.patch_info_list[idx]
        patch_data = self.image[
            patch_info[0] : patch_info[0] + self.patch_size[0],
            patch_info[1] : patch_info[1] + self.patch_size[1],
        ]
        if self.preproc is not None:
            patch_data = self.preproc(patch_data)
        return patch_data, patch_info
