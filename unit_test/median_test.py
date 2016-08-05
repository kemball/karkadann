from puremedian import *
from database import *
import unittest as ut
from time import time

class medianTest(ut.TestCase):


	def test_speed(self):
		seqa = "DDVLLVHGGTPLEGEIRVRGAKNLVPKAMVAALLGSGPSRLRNVPDIRDVRVVRGLLQLHGVTVRPGDEAGELILDPSHVESANVADIDAHAGSSRIPILFCGPLLHRLGHAFIPGLGGCDIGGRPVDFHFDVLRQFGATIEKRADGQYLEAPQRLRGCKIRLPYPSVGSTEQVLLTAVLAEGVTELSNAAVEPEIEDLICVLQKMGAIISMDTDRTIRITGVDKLGGYNHRALPDRLEAASWASAALATEGNIYVRGAQQRSMMTFLNTYRKVGGAFEIDDEGIRFWHPGGSLDAIALETDVHPGFQTDWQQPLVVALTQASGLSIVHETVYESRLGFTSALNQMGAHIQLYRECLGGSACRFGQRNFLHSAVVSGPTKLQGADLVIPDLRGGFSYLIAALAAQGTSRVHGIDLINRGYENFMEKLVSLGAHV"
		seqb= "WESKRIILAKERVGFSLHETILYAGTETSMWYANHIEAVLCVEGEAELTNDETGEKHTITPGTMYLLDGHEKHTMRIKEDFRCVCVFNPPVTGREDHDENGVYP"
		self.assertEqual(seq_identity(seqa,seqb),.14221218961625282)
		before = time()
		for i in range(0,10):
			seq_identity(seqa,seqb)
		print "seq_identity takes %s ms" %(time()-before)

	def test_runs(self):
		before = time()
		with get_cursor() as cur:
			cur.execute("select distinct id from clusters where classification !='cf_putative' limit 2;")
			clusters = list(Cluster.get_many([x for (x,) in cur.fetchall()]))
		self.assertTrue(clusters[0].is_real(),clusters[1].is_real())
		pure_median_identity(clusters[0],clusters[1])
		print "1 run takes %s seconds " %(time()-before)

	def test_reflexive(self):
		with get_cursor() as cur:
			cur.execute("select distinct id from clusters where classification !='cf_putative' limit 10;")
			clusters = list(Cluster.get_many([x for (x,) in cur.fetchall()]))
		forward = pure_median_identity(clusters[0],clusters[1])
		reverse = pure_median_identity(clusters[1],clusters[0])
		self.assertEqual(forward,reverse)