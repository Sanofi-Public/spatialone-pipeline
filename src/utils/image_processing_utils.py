import numpy as np


class ImageManipulation:
    """
    Class that can be utilized in any part of the pipeline for image processing
    """

    def __init__(self, img_arr, constant_value):
        self.img_shape = img_arr.shape
        self.constant_value = constant_value

    def crop_and_pad_img(
        self, img_arr: np.ndarray, img_shape: tuple, dim_length: int, dim_width: int
    ) -> np.ndarray:
        """Resize image.
        Arguments:
            img_arr {np.ndarray} -- source image
            img_shape {Tuple[int]} -- image shape
            dim_length {int} -- resized length dimension
            dim_width {int} -- resized length dimension
        Returns:
            np.ndarray -- result resized np.ndarray
        """
        img_length = img_shape[0]
        img_width = img_shape[1]
        if dim_width > img_width and dim_length < img_length:
            padded_img = self.padding(
                img_arr, img_length, dim_width, constant_value=self.constant_value
            )
            ldif = (img_length - dim_length) // 2
            if img_arr.ndim == 3:
                resized_img = padded_img[ldif : ldif + dim_length, :, :]
            else:
                resized_img = padded_img[ldif : ldif + dim_length, :]
        elif dim_length > img_length and dim_width < img_width:
            padded_img = self.padding(
                img_arr, dim_length, img_width, constant_value=self.constant_value
            )
            wdif = (img_width - dim_width) // 2
            if img_arr.ndim == 3:
                resized_img = padded_img[:, wdif : wdif + dim_width, :]
            else:
                resized_img = padded_img[:, wdif : wdif + dim_width]
        elif dim_length > img_length and dim_width > img_width:
            resized_img = self.padding(
                img_arr, dim_length, dim_width, constant_value=self.constant_value
            )
        else:
            padded_img = img_arr
            wdif = (img_width - dim_width) // 2
            ldif = (img_length - dim_length) // 2
            if padded_img.ndim == 3:
                resized_img = padded_img[
                    ldif : ldif + dim_length, wdif : wdif + dim_width, :
                ]
            else:
                resized_img = padded_img[
                    ldif : ldif + dim_length, wdif : wdif + dim_width
                ]
        return resized_img

    def padding(self, array: np.array, x_dim: int, y_dim: int, constant_value: int):
        """Pad around input array.
        Arguments:
            array {np.array} -- any numpy array
            x_dim {int} : desired height
            y_dim {int}: desirex width
            contant_value {int}: value to pad the array with
        Returns:
            np.array -- padded array
        """

        h = array.shape[0]
        w = array.shape[1]

        a = (x_dim - h) // 2
        aa = x_dim - a - h

        b = (y_dim - w) // 2
        bb = y_dim - b - w

        if array.ndim == 3:
            return np.pad(
                array,
                pad_width=((a, aa), (b, bb), (0, 0)),
                mode="constant",
                constant_values=constant_value,
            )
        else:
            return np.pad(
                array,
                pad_width=((a, aa), (b, bb)),
                mode="constant",
                constant_values=constant_value,
            )
