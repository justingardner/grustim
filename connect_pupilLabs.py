import nest_asyncio
from pupil_labs.realtime_api.simple import Device
import time

nest_asyncio.apply()
device = Device(address=ip, port= "8080")
