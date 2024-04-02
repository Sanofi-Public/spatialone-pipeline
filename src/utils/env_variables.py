"""
Env_variables for .env
"""

import os

from dotenv import load_dotenv


class Singleton:
    """This is a singleton class"""

    _instance = None

    def __new__(cls, test_mode: bool = False, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(Singleton, cls).__new__(
                cls,
                *args,
                **kwargs,
            )
        return cls._instance


class Env(Singleton):
    """Class definition to get the Environment Variables from .env"""

    def __init__(self, test_mode: bool = False):
        """Load the variables from .env file and store them
        in the private variables of Env class
        """
        load_dotenv(override=True)

        self._project_name = os.environ.get("PROJECT_NAME")

    # Creating functions inside the class to retrive the variables in another script
    # project specific functions
    @property
    def project_name(self):
        """x"""
        return self._project_name


if __name__ == "__main__":
    e = Env()
