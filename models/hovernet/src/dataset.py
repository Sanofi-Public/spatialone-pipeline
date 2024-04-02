""" dataset.py
"""
import glob
import json

import cv2
import numpy as np
import scipy.io as sio
import skimage.measure
import tifffile as tiff

NOT_SUPPORT_STRING = "Not support"


class __AbstractDataset(object):
    """Abstract class for interface of subsequent classes.
    Main idea is to encapsulate how each dataset should parse
    their images and annotations.

    """

    def load_img(self, path):
        """load_img

        Args:
            path (_type_): _description_

        Raises:
            NotImplementedError: _description_
        """
        raise NotImplementedError

    def load_ann(self, path, with_type=False):
        """_summary_

        Args:
            path (_type_): _description_
            with_type (bool, optional): _description_. Defaults to False.

        Raises:
            NotImplementedError: _description_
        """
        raise NotImplementedError


####
class __Kumar(__AbstractDataset):
    """Defines the Kumar dataset as originally introduced in:

    Kumar, Neeraj, Ruchika Verma, Sanuj Sharma, Surabhi Bhargava, Abhishek Vahadane,
    and Amit Sethi. "A dataset and a technique for generalized nuclear segmentation for
    computational pathology." IEEE transactions on medical imaging 36, no. 7 (2017): 1550-1560.

    """

    def load_img(self, path):
        """load_img

        Args:
            path (_type_): _description_

        Returns:
            _type_: _description_
        """
        return cv2.cvtColor(cv2.imread(path), cv2.COLOR_BGR2RGB)

    def load_ann(self, path, with_type=False):
        """load_ann

        Args:
            path (_type_): _description_
            with_type (_type_, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        # assumes that ann is HxW
        assert not with_type, NOT_SUPPORT_STRING
        ann_inst = sio.loadmat(path)["inst_map"]
        ann_inst = ann_inst.astype("int32")
        ann = np.expand_dims(ann_inst, -1)
        return ann


####
class __CPM17(__AbstractDataset):
    """Defines the CPM 2017 dataset as originally introduced in:

    Vu, Quoc Dang, Simon Graham, Tahsin Kurc, Minh Nguyen Nhat To, Muhammad Shaban,
    Talha Qaiser, Navid Alemi Koohbanani et al. "Methods for segmentation and classification
    of digital microscopy tissue images." Frontiers in bioengineering and biotechnology 7 (2019).

    """

    def load_img(self, path):
        """load_img

        Args:
            path (_type_): _description_

        Returns:
            _type_: _description_
        """
        return cv2.cvtColor(cv2.imread(path), cv2.COLOR_BGR2RGB)

    def load_ann(self, path, with_type=False):
        """load_ann

        Args:
            path (_type_): _description_
            with_type (_type_, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        assert not with_type, NOT_SUPPORT_STRING
        # assumes that ann is HxW
        ann_inst = sio.loadmat(path)["inst_map"]
        ann_inst = ann_inst.astype("int32")
        ann = np.expand_dims(ann_inst, -1)
        return ann
        assert not with_type, NOT_SUPPORT_STRING
        # assumes that ann is HxW
        ann_inst = sio.loadmat(path)["inst_map"]
        ann_inst = ann_inst.astype("int32")
        ann = np.expand_dims(ann_inst, -1)
        return ann


####
class __CoNSeP(__AbstractDataset):
    """Defines the CoNSeP dataset as originally introduced in:

    Graham, Simon, Quoc Dang Vu, Shan E. Ahmed Raza, Ayesha Azam, Yee Wah Tsang, Jin Tae Kwak,
    and Nasir Rajpoot. "Hover-Net: Simultaneous segmentation and classification of nuclei in
    multi-tissue histology images." Medical Image Analysis 58 (2019): 101563

    """

    def load_img(self, path):
        """load_img

        Args:
            path (_type_): _description_

        Returns:
            _type_: _description_
        """
        return cv2.cvtColor(cv2.imread(path), cv2.COLOR_BGR2RGB)

    def load_ann(self, path, with_type=False):
        """load_ann

        Args:
            path (_type_): _description_
            with_type (_type_, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        assert not with_type, NOT_SUPPORT_STRING
        # assumes that ann is HxW
        ann_inst = sio.loadmat(path)["inst_map"]
        # assumes that ann is HxW
        ann_inst = sio.loadmat(path)["inst_map"]
        if with_type:
            ann_type = sio.loadmat(path)["type_map"]

            # merge classes for CoNSeP (in paper we only utilise 3 nuclei classes and background)
            # If own dataset is used, then the below may need to be modified
            ann_type[(ann_type == 3) | (ann_type == 4)] = 3
            ann_type[(ann_type == 5) | (ann_type == 6) | (ann_type == 7)] = 4

            ann = np.dstack([ann_inst, ann_type])
            ann = ann.astype("int32")
        else:
            ann = np.expand_dims(ann_inst, -1)
            ann = ann.astype("int32")

        return ann


class __Lymph(__AbstractDataset):
    """
    Sample method to load Internal Sanofi Manually Annotated Data
    """

    def load_img(self, image_path):
        """
        Loading images in the form of tiffs
        """
        raw_img = cv2.cvtColor(tiff.imread(image_path), cv2.COLOR_BGR2RGB)
        return raw_img

    def load_manual_annotations(self, annotation_path):
        """
        Load in manual annotations from geojson file
        """
        with open(annotation_path) as annotations:
            file_contents = annotations.read()
            self.expert_annotations = json.loads(file_contents)
            if type(self.expert_annotations) == dict:
                self.expert_annotations = self.expert_annotations["features"]
        return self.expert_annotations

    def manual_annotation_polygons(self, expert_annotations):
        """
        Gather cell annotation coordinates from annotation file
        """
        self.cell_polygons = {}
        for item in self.expert_annotations:
            cell_type_polygons = []
            properties = item.get("properties")
            classification = properties.get("classification")
            label = classification.get("name")
            geometry = item.get("geometry")
            coordinates = geometry.get("coordinates")
            pts = np.array(coordinates[0], np.int32).reshape(-1, 1, 2)
            cell_type_polygons.append(pts)
            if label in self.cell_polygons.keys():
                self.cell_polygons[label].append(pts)
            else:
                self.cell_polygons[label] = cell_type_polygons
        return self.cell_polygons

    def rgb2label(self, img, color_codes=None, one_hot_encode=False):
        """
        Convert rgb image to label map
        """
        if color_codes is None:
            color_codes = {
                val: i
                for i, val in enumerate(set(tuple(v) for m2d in img for v in m2d))
            }
        n_labels = len(color_codes)
        result = np.ndarray(shape=img.shape[:2], dtype=int)
        result[:, :] = -1
        for rgb, idx in color_codes.items():
            result[(img == rgb).all(2)] = idx

        if one_hot_encode:
            one_hot_labels = np.zeros((img.shape[0], img.shape[1], n_labels))
            # one-hot encoding
            for c in range(n_labels):
                one_hot_labels[:, :, c] = (result == c).astype(int)
            result = one_hot_labels

        return result, color_codes

    def create_legend(self, type_info_path):
        """create legend from json file"""
        with open(type_info_path) as json_file:
            type_info = json.load(json_file)
        color_to_label = {}
        cell_type_legend = {}
        for key, value in type_info.items():
            cell_type_legend[value[0]] = tuple(value[1])
            color_to_label[tuple(value[1])] = key
        return color_to_label, cell_type_legend

    def draw_annotations(
        self, cell_polygons, raw_img, color_to_label, cell_type_legend, with_type=True
    ):
        """
        Draw annotations from geojson files onto raw image and create instance and type maps
        """
        ##add way to load in legend from json
        ann_mask = raw_img.copy()
        for key, value in cell_polygons.items():
            if key in cell_type_legend.keys():
                color = cell_type_legend.get(key)
            else:
                color = cell_type_legend.get("nolabel")
            cv2.fillPoly(ann_mask, value, color=color)
        type_mask, _ = self.rgb2label(ann_mask, color_to_label)
        inst_mask = skimage.measure.label(type_mask)
        type_mask[type_mask > 10] = 10
        if with_type:
            ann = np.dstack([inst_mask, type_mask])
            ann = ann.astype("int32")
        else:
            ann = np.expand_dims(inst_mask, -1)
            ann = ann.astype("int32")
        return ann

    def load_ann(
        self, annotation_path=None, type_info_path=None, raw_img=None, with_type=None
    ):
        """
        Load annotation data
        """
        ##add way to load in legend from json
        expert_annotations = self.load_manual_annotations(annotation_path)
        cell_polygons = self.manual_annotation_polygons(expert_annotations)
        color_to_label, cell_type_legend = self.create_legend(type_info_path)
        ann = self.draw_annotations(
            cell_polygons,
            raw_img,
            color_to_label,
            cell_type_legend,
            with_type=with_type,
        )
        return ann


####
def get_dataset(name):
    """Return a pre-defined dataset object associated with `name`."""
    name_dict = {
        "kumar": lambda: __Kumar(),
        "cpm17": lambda: __CPM17(),
        "consep": lambda: __CoNSeP(),
        "lymphocytes": lambda: __Lymph(),
    }
    if name.lower() in name_dict:
        return name_dict[name]()
    else:
        assert False, "Unknown dataset `%s`" % name
