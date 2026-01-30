import logging

# Create a logger for the package
logger = logging.getLogger("icepop")

# Set the log level (e.g., DEBUG, INFO, WARNING, ERROR, CRITICAL)
logger.setLevel(logging.DEBUG)

# Create a console handler
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)

# Create a formatter and set it for the handler
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console_handler.setFormatter(formatter)

# Add the handler to the logger
logger.addHandler(console_handler)

# Prevent propagation to the root logger (optional)
logger.propagate = False
