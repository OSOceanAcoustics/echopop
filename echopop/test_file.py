class CONFIG_NEW(InputModel, arbitrary_types_allowed=True):
    test: str

    def __init__(self, filename, **kwargs):
        # self.model_config["title"] = filename
        # self.__class__.__name__ = filename

        try:
            super().__init__(**kwargs)
        except ValidationError as e:
            new_message = str(e).replace(
                self.__class__.__name__, 
                f"configuration parameters defined in {filename}"
            )
            raise ValueError(new_message) from e

        # try:
        #     super().__init__(**kwargs)
        # except ValidationError as e:
        #     raise e

        
CONFIG_NEW(filename="t", **dict(test="a")).model_config
CONFIG_NEW(filename="t", **dict()).__str__()
CONFIG_NEW(filename="t", **dict())