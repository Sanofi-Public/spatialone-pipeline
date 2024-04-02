"""
Pre-processing functions for image_seg
"""
import cv2
import numpy as np
import pandas as pd
import skimage

# import staintools
from patchify import patchify, unpatchify
from skimage import morphology
from skimage.measure import block_reduce

from src.utils.image_processing_utils import ImageManipulation


class Preprocessing:
    """
    Class to preprocess input H&E image
    """

    def __init__(self, img_arr, patch_size, overlap):
        self.img_shape = img_arr.shape
        self.patch_size = patch_size
        self.overlap = overlap
        self.downsampled_shape = None

    def resize_img_dim(self, downsampled_shape: tuple) -> int:
        """Find dimension to resize to for patching to work properly.
        Arguments:
            downsampled_shape {Tuple[int]} -- downsampled image shape
            patch_size {int} -- patch size
            overlap {int} -- step size
        Returns:
            int -- dimension to resize image to
        """
        maximum_width = max(
            downsampled_shape[0], downsampled_shape[1]
        )  # if largest dimension is taken, smaller side will need padding
        max_bound = maximum_width // self.overlap
        dim = ((max_bound - 1) * self.overlap) + self.patch_size
        print("Resized dimension: ", dim)
        return dim

    def resize_img(
        self, img_arr: np.ndarray, dim_length: int, dim_width: int
    ) -> np.ndarray:
        """Resize image.
        Arguments:
            img_arr {np.ndarray} -- source image
            dim {int} -- resized dimension
        Returns:
            np.ndarray -- result resized np.ndarray
        """
        image_manipulator = ImageManipulation(img_arr, 255)
        resized_img = image_manipulator.crop_and_pad_img(
            img_arr, image_manipulator.img_shape, dim_length, dim_width
        )
        return resized_img

    def detect_tissue_bbox(self, pos_csv: pd.DataFrame) -> list:
        """Resize image.
        Arguments:
            pos_csv{pd.dataframe} -- tissue positions list table
        Returns:
            list -- list of bbox coordinates [x1,y1,x2,y2]
        """
        pos_csv.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ]
        pos_csv = pos_csv[pos_csv.in_tissue != 0]
        pos_csv = pos_csv.drop(pos_csv.index[pos_csv["pxl_col_in_fullres"] < 0])
        pos_csv = pos_csv.drop(pos_csv.index[pos_csv["pxl_row_in_fullres"] < 0])
        pos_csv = pos_csv.drop(
            pos_csv.index[pos_csv["pxl_row_in_fullres"] > self.img_shape[0]]
        )
        pos_csv = pos_csv.drop(
            pos_csv.index[pos_csv["pxl_col_in_fullres"] > self.img_shape[1]]
        )
        max_width = pos_csv.max(axis=0)["pxl_col_in_fullres"] + self.patch_size
        max_length = pos_csv.max(axis=0)["pxl_row_in_fullres"] + self.patch_size
        min_width = pos_csv.min(axis=0)["pxl_col_in_fullres"] - self.patch_size
        min_length = pos_csv.min(axis=0)["pxl_row_in_fullres"] - self.patch_size

        x0 = max(0, min_width)
        x1 = min(self.img_shape[1], max_width)
        y0 = max(0, min_length)
        y1 = min(self.img_shape[0], max_length)
        bbox = [y0, y1, x0, x1]
        return bbox

    def pad_all_around(self, img_arr: np.ndarray):
        """Pads the image with patch_size pixels all around

        Args:
            img_arr (np.ndarray): _description_
        Returns:
            _type_: _description_
        """
        img_arr = np.pad(
            img_arr,
            (
                (self.patch_size, self.patch_size),
                (self.patch_size, self.patch_size),
                (0, 0),
            ),
            "constant",
            constant_values=(255,),
        )
        return img_arr

    def undo_padding(self, mask):
        """Undo the padding that was added

        Args:
            mask (_type_): _description_

        Returns:
            _type_: _description_
        """
        mask = mask[
            self.patch_size : -self.patch_size, self.patch_size : -self.patch_size
        ]
        return mask

    def adjust_tissue_boundary(self, boundary):
        """Shifts the tissue boundaries ([y0, y1, x0, x1]) to
        account for the padding added

        Args:
            boundary (_type_): _description_

        Returns:
            _type_: _description_
        """
        boudary_adj = [i + self.patch_size for i in boundary]
        return boudary_adj

    def crop_tissue_region(
        self, img_arr: np.ndarray, bbox: list[int], patch_size: int, pad: bool = True
    ) -> np.ndarray:
        """Crop Image based on tissue detected coordinates image.
        Arguments:
            img_arr {np.ndarray} -- source image
            bbox {List[int]} -- bounding box coordinates of tissue detected region
            patch_size {int} -- patch size
        Returns:
            np.ndarray -- result cropped np.ndarray
        """
        img = img_arr[
            bbox[0] - (patch_size) : bbox[1] + (patch_size),
            bbox[2] - (patch_size) : bbox[3] + (patch_size),
        ]
        return img

    def tissue_detector(self, img_arr: np.ndarray) -> np.ndarray:
        """
        Isolate Tissue from the image and remove artifacts in the background.
        Arguments:
            img_arr {np.ndarray} -- source image
        Returns:
            np.ndarray -- result cropped np.ndarray
        """
        # Convert to Grayscale
        gray = cv2.cvtColor(img_arr, cv2.COLOR_RGB2GRAY)
        # Otsu Threshold to remove background artifact
        ret, thresh1 = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        mask = morphology.remove_small_objects(
            thresh1 == 0, min_size=16 * 16, connectivity=2
        )
        mask = morphology.remove_small_holes(mask, area_threshold=128 * 128)
        mask = morphology.binary_dilation(mask, morphology.disk(16))
        # Make background empty on original image
        img_arr[np.where(mask == 0)] = 255
        return img_arr

    def norm_stain_remove(self, img_arr: np.ndarray) -> np.ndarray:
        """Normalize H&E image.
        Arguments:
            img_arr {np.ndarray} -- source image
        Returns:
            np.ndarray -- normalized np.ndarray
        """
        # Macenko removed due to bug causing entire image to become empty
        hed_img = skimage.color.rgb2hed(img_arr)
        # img = hed_img.astype("uint8")
        # normalizer = staintools.StainNormalizer(method="Macenko")
        # normalizer.fit(img)
        # transformed = normalizer.transform(img)
        return hed_img

    def downsample(self, img_arr: np.ndarray, downsample_factor: int) -> np.ndarray:
        """Downsample image
        Arguments:
            img_arr {np.ndarray} -- source image
            downsample_factor {int} -- downsampling factor
        Returns:
            np.ndarray -- result downsampled np.ndarray
        """
        r = block_reduce(
            img_arr[:, :, 0], (downsample_factor, downsample_factor), np.mean
        )
        g = block_reduce(
            img_arr[:, :, 1], (downsample_factor, downsample_factor), np.mean
        )
        b = block_reduce(
            img_arr[:, :, 2], (downsample_factor, downsample_factor), np.mean
        )
        downsampled_img = np.stack((r, g, b), axis=-1)
        return downsampled_img

    def generate_patches(self, img_arr: np.ndarray, n_channels: int) -> np.ndarray:
        """Split image into pathches.
        Arguments:
            img_arr {np.ndarray} -- source image
            patch_size {int} -- patch size
            n_channels {int} -- number of channels of images
            overlap {int} -- overlap
        Raises:
            NotImplementedError: Raises error while if the image shape is not supported
        Returns:
            np.ndarray -- result np.ndarray with patches
        """
        if len(img_arr.shape) == 3:
            return patchify(
                img_arr, (self.patch_size, self.patch_size, n_channels), self.overlap
            )
        elif len(img_arr.shape) == 2:
            return patchify(img_arr, (self.patch_size, self.patch_size), self.overlap)
        else:
            raise NotImplementedError(
                "Only 2D images or 3D (2D RGB) images are supported."
            )
