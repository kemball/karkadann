-- MySQL dump 10.13  Distrib 5.5.47, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: karkadann
-- ------------------------------------------------------
-- Server version	5.5.47-0ubuntu0.14.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `assemblies`
--

DROP TABLE IF EXISTS `assemblies`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `assemblies` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `gb_record` varchar(255) DEFAULT NULL,
  `fastq` varchar(255) DEFAULT NULL,
  `assembled` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `genome_id` int(11) NOT NULL,
  `accession` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `unique_file` (`gb_record`),
  KEY `genome_id` (`genome_id`),
  CONSTRAINT `assembly_genome` FOREIGN KEY (`genome_id`) REFERENCES `genomes` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=266 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `assemblies`
--

LOCK TABLES `assemblies` WRITE;
/*!40000 ALTER TABLE `assemblies` DISABLE KEYS */;
/*!40000 ALTER TABLE `assemblies` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `contigs`
--

DROP TABLE IF EXISTS `contigs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `contigs` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `sequence` longtext,
  `assembly_id` int(11) NOT NULL,
  `accession` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `assembly_id` (`assembly_id`),
  CONSTRAINT `contig_assembly` FOREIGN KEY (`assembly_id`) REFERENCES `assemblies` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=2689 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `contigs`
--

LOCK TABLES `contigs` WRITE;
/*!40000 ALTER TABLE `contigs` DISABLE KEYS */;
/*!40000 ALTER TABLE `contigs` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `genes`
--

DROP TABLE IF EXISTS `genes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genes` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `translation` text,
  `start` varchar(100) DEFAULT NULL,
  `end` varchar(100) DEFAULT NULL,
  `strand` enum('1','-1') DEFAULT NULL,
  `contig` int(11) NOT NULL,
  `accession` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `contig` (`contig`),
  CONSTRAINT `genes_contig` FOREIGN KEY (`contig`) REFERENCES `contigs` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=124251 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `genes`
--

LOCK TABLES `genes` WRITE;
/*!40000 ALTER TABLE `genes` DISABLE KEYS */;
/*!40000 ALTER TABLE `genes` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `genomes`
--

DROP TABLE IF EXISTS `genomes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genomes` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `added` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `unique_name` (`name`)
) ENGINE=InnoDB AUTO_INCREMENT=77228 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `genomes`
--

LOCK TABLES `genomes` WRITE;
/*!40000 ALTER TABLE `genomes` DISABLE KEYS */;
/*!40000 ALTER TABLE `genomes` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `genus_species`
--

DROP TABLE IF EXISTS `genus_species`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genus_species` (
  `genome_id` int(11) DEFAULT NULL,
  `binomial` text,
  KEY `genome_id` (`genome_id`),
  CONSTRAINT `genus_species_ibfk_1` FOREIGN KEY (`genome_id`) REFERENCES `genomes` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `genus_species`
--

LOCK TABLES `genus_species` WRITE;
/*!40000 ALTER TABLE `genus_species` DISABLE KEYS */;
/*!40000 ALTER TABLE `genus_species` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `hits`
--

DROP TABLE IF EXISTS `hits`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hits` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `score` float DEFAULT NULL,
  `hmm` varchar(100) DEFAULT NULL,
  `gene` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `unique_hmm` (`hmm`,`gene`),
  KEY `gene_map` (`gene`),
  CONSTRAINT `gene_map` FOREIGN KEY (`gene`) REFERENCES `genes` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `hits`
--

LOCK TABLES `hits` WRITE;
/*!40000 ALTER TABLE `hits` DISABLE KEYS */;
/*!40000 ALTER TABLE `hits` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `hmm_storage`
--

DROP TABLE IF EXISTS `hmm_storage`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hmm_storage` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `filename` varchar(200) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `hmm_storage`
--

LOCK TABLES `hmm_storage` WRITE;
/*!40000 ALTER TABLE `hmm_storage` DISABLE KEYS */;
/*!40000 ALTER TABLE `hmm_storage` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `names`
--

DROP TABLE IF EXISTS `names`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `names` (
  `genome_id` int(11) NOT NULL,
  `cc` varchar(50) DEFAULT NULL,
  `name` varchar(255) DEFAULT NULL,
  KEY `genome_id` (`genome_id`),
  CONSTRAINT `names_ibfk_1` FOREIGN KEY (`genome_id`) REFERENCES `genomes` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `names`
--

LOCK TABLES `names` WRITE;
/*!40000 ALTER TABLE `names` DISABLE KEYS */;
/*!40000 ALTER TABLE `names` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2016-03-18 15:36:34
