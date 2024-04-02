"""
Post-processing functions for image_seg
"""
import itertools
import math
import warnings
from typing import Dict, List, Tuple

import cv2
import numpy as np
from skimage.measure import label

from src.utils.image_processing_utils import ImageManipulation


def stitching_instance_segmentation(
    patches: np.array,
    patch_size: Tuple[int],
    step_size: int,
    whole_size: Tuple[int],
    minimum_pixels: int = 4,
) -> np.array:
    """This function is adapted from
    https://github.com/michaellee1/CellSeg/blob/master/src/cvstitch.py
    Lee, M.Y., Bedia, J.S., Bhate, S.S. et al. CellSeg: a robust, pre-trained
    nucleus segmentation and pixel quantification software for highly multiplexed
    fluorescence images. BMC Bioinformatics 23, 46 (2022).
    https://doi.org/10.1186/s12859-022-04570-9
    Some notations:
    - A mask is an instance segmentation mask that may contains more than 1
     cells segmentation results.
    - An instance refers to the region (usually a connected composant) on
     an instance segmentation mask, whose pixels share same pixel value,
     i.e. the instance id.
    - A patch is a crop of whole image/mask/slide, in this funciton, an instance
     segmentation mask crop (for a image crop).
    Args:
        patches: patches of shape (nrows, ncols, patch_h, patch_w)
        patch_size: the dimension of single patch (patch_h, patch_w)
        step_size: overlap between patches in terms of number of pixels.
        whole_size: the whole image/mask size.
        minimum_pixels (int): minimum number of pixels for one instance
        to not be considered as artifacts.
    Returns:
        numpy.ndarray: whole mask stitched based on input patches
    """
    # I. prepare the instance segmetation masks with right format and indexing
    patch_h, patch_w = patch_size

    nrows = patches.shape[0]
    ncols = patches.shape[1]
    # check inputs
    if (
        (whole_size[0] - patch_h) // step_size != nrows - 1
        or (whole_size[1] - patch_w) // step_size != ncols - 1
        or patches.shape[-2] != patch_size[0]
        or patches.shape[-1] != patch_size[1]
    ):
        raise ValueError(
            "The shape of input patches does not match patching and size parameters."
        )

    # Assign unique instances IDs across the full image:
    # If an instance appears in multiple tiles, reindexing would assign
    # this instance with multiple IDs, which, in the following steps, would be
    # populated into different layers and identified as overlapped instances.

    patches = reindex_masks(patches)

    # II. remove cells under a certain size to get rid of artifacts
    print(f"Removing instances/cells with area less than {minimum_pixels} pixels.")
    instance_ids = []
    instance_sizes = []
    for idrow, idcol in itertools.product(range(nrows), range(ncols)):
        patches[idrow, idcol], sizes, ids = remove_small_instances(
            patches[idrow, idcol], threshold=minimum_pixels
        )
        instance_sizes.extend(sizes)
        instance_ids.extend(ids)

    size_dict = dict(zip(instance_ids, instance_sizes))
    del instance_ids, instance_sizes

    # III. populate an expanded mask of same dimension as whole image/mask/slide
    # with multiple layers. Currently 4 layers are used, which ensure good
    # stitching for patches with <= 50% overlapping.

    nlayers = 2
    expanded_mask_arr = populate_whole_mask_with_patches(
        patches, whole_size, step_size, nlayers=nlayers
    )

    # IX. maximum overlapping and choose the biggest instances if duplicated,
    # return a stitched output of shape (whole_h, whole_w) combining all layers.

    # Deal with overlapping: get conflicted instances and compare mask sizes
    # of overlapping masks, only keeping the largest mask in each conflict and
    # removing all other masks.
    expanded_mask_arr = resolve_conflict(expanded_mask_arr, size_dict=size_dict)

    # relabel in case of indexing issues
    expanded_mask_arr = label(expanded_mask_arr)

    return expanded_mask_arr


def reindex_single_mask(mask: np.array) -> np.array:
    """Reindex a series of instance segmentation masks
    If the total number of instances equals to the maximum
    instance index, we consider that the instance segmentation
    mask is well indexed and return a copy of input mask.
    Args:
        mask (numpy.ndarray): single instance segmentation mask
    Returns:
        numpy.ndarray: instance segmentation mask with all instances
         numbered from 1, 2, .... to num_instances found on input mask.
    """
    # reindex instance on current mask to number them from 1, 2...
    if len(np.unique(mask)) != mask.max() + 1 or (mask < 0).any():
        unique_indices = np.sort(np.unique(mask))
        new_indices = np.arange(len(unique_indices))
        mapping = np.digitize(mask.ravel(), unique_indices, right=True)
        reindex_mask = new_indices[mapping].reshape(mask.shape)
        return reindex_mask
    else:
        return mask


def reindex_masks(masks: np.array) -> np.array:
    """Reindex a series of instance segmentation masks
    Args:
        masks (numpy.ndarray): input instance segmentations must be
         of shape (nrow, ncol, height, width). The instances on a
         single instance segmentation mask have unique index (not necessarily
         starting from 1) but all instances from all masks may share
         duplicated index.
    Returns:
        numpy.ndarray: all instances have unique index among
         all masks from 1, 2, .... to num_instances.
    """
    current_index = 0
    reindexed_masks = masks.copy()

    for idrow, idcol in itertools.product(range(masks.shape[0]), range(masks.shape[1])):
        mask = masks[idrow, idcol]
        if np.any(mask.astype(bool)):
            # reindex instances on current mask from 1, 2... if needed
            reindex_mask = reindex_single_mask(mask)

            # add current_index / count of instances from previous masks
            reindex_mask[mask > 0] += current_index
            current_index += len(np.unique(reindex_mask)) - 1

            # update output
            reindexed_masks[idrow, idcol] = reindex_mask
    return reindexed_masks


def remove_small_instances(
    mask: np.array, threshold: float = 8
) -> Tuple[np.array, List[int], List[int]]:
    """This function remove instances with size (in pixels) smaller
    than a threshold from instance segmentation results.
    Args:
        mask (numpy.ndarray): instance segmentation mask or masks
        threshold (float): minimum number of pixels for one instance to not
        be considered as artifacts.
    Returns:
        numpy.ndarray: filtered instance segmentation mask or masks
        List[int]: List of pixel sizes for kept instances
        List[int]: List of instance ids for kept instances (note that
         instance 0 / background can be included)
    """
    filtered_mask = mask.copy()
    instance_ids, sizes = np.unique(mask, return_counts=True)
    keep_indices = list(sizes >= threshold)
    for instance_id, keep in zip(instance_ids, keep_indices):
        if not keep:
            filtered_mask[mask == instance_id] = 0

    return filtered_mask, list(sizes[keep_indices]), list(instance_ids[keep_indices])


def populate_whole_mask_with_patches(
    patches: np.array, whole_size: Tuple[int], step_size: int, nlayers: int = 4
) -> np.array:
    """This function populates an expanded mask of same dimension as whole image/mask/slide
    with multiple layers
    Args:
        patches: patches of shape (nrows, ncols, patch_h, patch_w)
        whole_size: the whole image/mask size.
        step_size: overlap between patches in terms of number of pixels.
        nlayers: layer of output expanded mask
    Returns:
        np.array: expanded mask of shape (whole_h, whole_w, nlayers)
         populated with pixel values from patches. An intermediate step
         for stitching.
    """

    # check inputs
    patch_h, patch_w = patches.shape[-2:]
    whole_h, whole_w = whole_size

    nrows = patches.shape[0]
    ncols = patches.shape[1]

    if (whole_size[0] - patch_h) // step_size != nrows - 1 or (
        whole_size[1] - patch_w
    ) // step_size != ncols - 1:
        raise ValueError("The shape of input patches and size parameters do not match.")

    level_overlapping = math.floor(math.sqrt(nlayers))
    if (
        step_size / patch_w < 1 / level_overlapping
        or step_size / patch_h < 1 / level_overlapping
    ):
        warnings.warn("Highly overlapping patches - recommend using a higher nlayers.")

    # initialize an empty expanded mask array

    expanded_mask_arr = np.zeros((whole_h, whole_w, nlayers), dtype=np.int32)

    # populate the expanded mask array with multiple layers
    layer_to_populate = 0

    curr_left = 0
    curr_top = 0
    for idrow in range(nrows):
        curr_left = 0
        for idcol in range(ncols):
            patch = patches[idrow, idcol]
            # make sure neighboring patches not populated into same layer
            layer_to_populate = (
                (idrow % level_overlapping) * level_overlapping + idcol
            ) % nlayers
            expanded_mask_arr[
                curr_top : (curr_top + patch_h),
                curr_left : (curr_left + patch_w),
                layer_to_populate,
            ] = patch
            curr_left += step_size
        curr_top += step_size

    return expanded_mask_arr


def check_conflict(
    expanded_mask_arr: np.array,
) -> np.array:
    """This function returns conflicted instances based on multi-layer populated
    mask.
    Args:
        expanded_mask_arr: expanded mask of shape (patch_h, patch_w, nlayers)
         populated with pixel values from patches. An intermediate step
         for stitching.
    Returns:
        np.array: summary of conflited instance ids arranged in an array of shape
         (nlayers, number_conflicts_found). For each conflict, the ids of instances
         overlapping found from all layers are recorded as a column in output array,
         which is filled with 0 if not found for some of the layers.
         Background (instance id=0) is excluded from the consideration of overlapping.
    """
    instance_overlaps = np.sum(expanded_mask_arr > 0, axis=2) > 1
    instance_overlaps_compress = np.zeros(
        (4, np.sum(instance_overlaps)), dtype=np.int64
    )

    for i in range(expanded_mask_arr.shape[2]):
        instance_overlaps_compress[i, :] = expanded_mask_arr[:, :, i][instance_overlaps]

    instance_conflicts = np.unique(instance_overlaps_compress, axis=1)
    return instance_conflicts


def resolve_conflict(
    expanded_mask_arr: np.array,
    size_dict: Dict[int, int],
) -> np.array:
    """This function returns stitched whole mask of shape (whole_h, whole_w)
    based on multi-layer populated mask of shape (whole_h, whole_w, nlayers).
    Args:
        expanded_mask_arr: expanded mask of shape (whole_h, whole_w, nlayers)
         populated with pixel values from patches. An intermediate step
         for stitching.
        size_dict: dictionary mapping instance ids to instance size in terms
         of number of pixels.
    Returns:
        np.array: stitched whole mask of shape (whole_h, whole_w)
    """
    # initialize a list for the ids of instances to remove (overlapping with other
    # instances but not the largest).
    instances_to_rem = []

    instance_conflicts = check_conflict(expanded_mask_arr)
    for i in range(instance_conflicts.shape[1]):
        overlap_instance_ids = instance_conflicts[:, i]
        # For example, `overlap_instance_ids`` can be (0, 66, 87, 0),
        # if number of layers is 4, which means we found instances 66 and 87
        # overlapping.
        # Note that 0 does not mean background but overlapping instances
        # not found on layer 0 and layer 3.
        overlap_instance_ids = [idi for idi in overlap_instance_ids if idi > 0]

        remove_indices = get_remove_indices_max_overlapping(
            overlap_instance_ids, size_dict
        )
        instances_to_rem.extend(remove_indices)

    # remove duplicated instances
    instances_to_rem = list(
        set(instances_to_rem)
    )  # get unique ids for instances to remove
    masklocs = np.isin(expanded_mask_arr, instances_to_rem)
    expanded_mask_arr[masklocs] = 0

    # sum up all layers to get an array of shape (whole_h, whole_w)
    full_mask_arr = np.sum(expanded_mask_arr, axis=2)

    # renumber the indices before returning array
    reindexed_full_arr = reindex_single_mask(full_mask_arr)

    # return stitched whole mask
    return reindexed_full_arr


def get_remove_indices_max_overlapping(
    instance_ids: np.array,
    size_dict: Dict[int, int],
) -> List[int]:
    """
    This function returns list of instances ids that are not the largest among
     overlapping indices.
    Args:
        instance_ids: an 1d array of instance ids.
        size_dict: dictionary mapping instance ids to instance size in terms
         of number of pixels.
    Returns:
        List[int]: list of instances ids that are not the largest among
         overlapping indices in `instance_ids`.
    """
    instance_sizes = [size_dict.get(inst_id) for inst_id in instance_ids]
    largest_instance_idx = np.argmax(np.array(instance_sizes))
    largest_instance = instance_ids[largest_instance_idx]
    remove_indices = [
        inst_id for inst_id in instance_ids if inst_id != largest_instance
    ]

    return remove_indices


def match_input_shape(
    segmentation: np.ndarray,
    input_im_shape: tuple,
    bbox: list,
    cropped_region: tuple,
    downsampled_factor: int,
    upsampled_only_flag=False,
) -> np.ndarray:
    """Reshape image to original shape
    Arguments:
        segmentation {np.ndarray} -- segmented mask
        input_im_shape {Tuple[int]} -- input image shape
        bbox {List[int]} -- bounding box coordinates of tissue detected region
        downsampled_factor {int} -- downsampled factor amount
    Returns:
        np.ndarray -- result input image size np.ndarray
    """
    upsample_factor = downsampled_factor
    if upsampled_only_flag is True:
        # Upsample mask
        upsample = cv2.resize(
            segmentation,
            None,
            fx=upsample_factor,
            fy=upsample_factor,
            interpolation=cv2.INTER_NEAREST,
        )
        print("upsample: ", upsample.shape)
        tissue_detected_length, tissue_detected_width = (bbox[1] - bbox[0]), (
            bbox[3] - bbox[2]
        )
        image_manipulator = ImageManipulation(upsample, 0)
        upsample_pad = image_manipulator.crop_and_pad_img(
            upsample,
            image_manipulator.img_shape,
            tissue_detected_length,
            tissue_detected_width,
        )
        return upsample_pad
    else:
        # Upsample mask
        upsample = cv2.resize(
            segmentation,
            None,
            fx=upsample_factor,
            fy=upsample_factor,
            interpolation=cv2.INTER_NEAREST,
        )
        print("upsample: ", upsample.shape)
        image_manipulator = ImageManipulation(upsample, 0)
        tissue_detected_length, tissue_detected_width = (bbox[1] - bbox[0]), (
            bbox[3] - bbox[2]
        )
        if (
            tissue_detected_length > input_im_shape[0]
            or tissue_detected_width > input_im_shape[1]
        ):
            print("upsampled image is larger than original cropped region")
            upsample_pad = image_manipulator.crop_and_pad_img(
                upsample,
                image_manipulator.img_shape,
                cropped_region[0],
                cropped_region[1],
            )
        else:
            upsample_pad = image_manipulator.crop_and_pad_img(
                upsample,
                image_manipulator.img_shape,
                tissue_detected_length,
                tissue_detected_width,
            )
        # check to see if size matches pre-downsampled size, otherwise crop
        print("upsample check shape: ", upsample_pad.shape)
        final_mask = np.zeros((input_im_shape[0], input_im_shape[1]))
        # Add to final_mask based on bbox coordinates
        final_mask[bbox[0] : bbox[1], bbox[2] : bbox[3]] = upsample_pad
    return final_mask, upsample_pad
