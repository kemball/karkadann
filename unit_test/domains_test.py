from karkadann.domains import *
from karkadann.domains import _call_mcl

from karkadann.database import *

import unittest as ut

class mclTest(ut.TestCase):
	def test_call_mcl(self):
		with get_cursor() as cur:
			cur.execute("select id from genes limit 1000;")
			records = [Gene(db_id=x).record for (x,) in cur.fetchall()]
			_call_mcl(records)