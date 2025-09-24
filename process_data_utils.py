import yaml

class ConfigLoader:
    """
    Clase para cargar y acceder a parámetros de configuración desde un archivo YAML.

    Uso:
        cfg = ConfigLoader("config.yaml")
        print(cfg.filter)  # Accede al parámetro 'filter'
        print(cfg.method)  # Accede al parámetro 'method' dentro de 'alignment'
    """

    def __init__(self, filepath="config.yaml"):
        """
        Inicializa la clase cargando el archivo YAML.

        Parámetros:
            filepath (str): Ruta al archivo YAML de configuración.
        """
        with open(filepath, "r") as file:
            self.config = yaml.safe_load(file)

    def get(self, key, default=None):
        """
        Accede a una clave de primer nivel del YAML.

        Parámetros:
            key (str): Clave a buscar.
            default: Valor por defecto si no se encuentra la clave.

        Retorna:
            Valor asociado a la clave o el valor por defecto.
        """
        return self.config.get(key, default)

    def __getattr__(self, name):
        """
        Permite acceder a cualquier parámetro del YAML como atributo.

        Si el parámetro está anidado (por ejemplo, dentro de 'alignment'),
        también se puede acceder directamente por su nombre.

        Ejemplo:
            cfg.method → accede a config['alignment']['method']
        """
        if name in self.config:
            return self.config[name]
        for key, value in self.config.items():
            if isinstance(value, dict) and name in value:
                return value[name]
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")