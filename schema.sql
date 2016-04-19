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
) ENGINE=InnoDB AUTO_INCREMENT=2054 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `assemblies`
--

LOCK TABLES `assemblies` WRITE;
/*!40000 ALTER TABLE `assemblies` DISABLE KEYS */;
/*!40000 ALTER TABLE `assemblies` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `cluster_names`
--

DROP TABLE IF EXISTS `cluster_names`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `cluster_names` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  `genome` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `genome` (`genome`),
  CONSTRAINT `cluster_names_ibfk_1` FOREIGN KEY (`genome`) REFERENCES `genomes` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=613 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `cluster_names`
--

LOCK TABLES `cluster_names` WRITE;
/*!40000 ALTER TABLE `cluster_names` DISABLE KEYS */;
/*!40000 ALTER TABLE `cluster_names` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `clusters`
--

DROP TABLE IF EXISTS `clusters`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `clusters` (
  `classification` varchar(100) NOT NULL,
  `gene` int(11) DEFAULT NULL,
  `id` int(11) DEFAULT NULL,
  UNIQUE KEY `unique_class` (`gene`,`id`),
  KEY `fk` (`id`),
  CONSTRAINT `clusters_ibfk_1` FOREIGN KEY (`gene`) REFERENCES `genes` (`id`) ON DELETE CASCADE,
  CONSTRAINT `clusters_ibfk_2` FOREIGN KEY (`id`) REFERENCES `cluster_names` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `clusters`
--

LOCK TABLES `clusters` WRITE;
/*!40000 ALTER TABLE `clusters` DISABLE KEYS */;
/*!40000 ALTER TABLE `clusters` ENABLE KEYS */;
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
) ENGINE=InnoDB AUTO_INCREMENT=44613 DEFAULT CHARSET=latin1;
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
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  `strand` int(11) DEFAULT NULL,
  `contig` int(11) NOT NULL,
  `accession` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `contig` (`contig`),
  CONSTRAINT `genes_contig` FOREIGN KEY (`contig`) REFERENCES `contigs` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=4484766 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=80659 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=237861 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `hits`
--

LOCK TABLES `hits` WRITE;
/*!40000 ALTER TABLE `hits` DISABLE KEYS */;
/*!40000 ALTER TABLE `hits` ENABLE KEYS */;
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

--
-- Table structure for table `orthogroups`
--

DROP TABLE IF EXISTS `orthogroups`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `orthogroups` (
  `batch` int(11) NOT NULL,
  `group_name` varchar(100) DEFAULT NULL,
  `gene` int(11) NOT NULL,
  UNIQUE KEY `bgene` (`batch`,`gene`),
  KEY `fkgene` (`gene`),
  CONSTRAINT `fkbatch` FOREIGN KEY (`batch`) REFERENCES `orthomcl_batches` (`id`) ON DELETE CASCADE,
  CONSTRAINT `fkgene` FOREIGN KEY (`gene`) REFERENCES `genes` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `orthogroups`
--

LOCK TABLES `orthogroups` WRITE;
/*!40000 ALTER TABLE `orthogroups` DISABLE KEYS */;
/*!40000 ALTER TABLE `orthogroups` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `orthomcl_batches`
--

DROP TABLE IF EXISTS `orthomcl_batches`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `orthomcl_batches` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `done` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=61 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `orthomcl_batches`
--

LOCK TABLES `orthomcl_batches` WRITE;
/*!40000 ALTER TABLE `orthomcl_batches` DISABLE KEYS */;
/*!40000 ALTER TABLE `orthomcl_batches` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `promer`
--

DROP TABLE IF EXISTS `promer`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `promer` (
  `score` float DEFAULT NULL,
  `l` int(11) DEFAULT NULL,
  `r` int(11) DEFAULT NULL,
  UNIQUE KEY `lr` (`l`,`r`),
  KEY `fkr` (`r`),
  CONSTRAINT `fkl` FOREIGN KEY (`l`) REFERENCES `cluster_names` (`id`) ON DELETE CASCADE,
  CONSTRAINT `fkr` FOREIGN KEY (`r`) REFERENCES `cluster_names` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `promer`
--

LOCK TABLES `promer` WRITE;
/*!40000 ALTER TABLE `promer` DISABLE KEYS */;
/*!40000 ALTER TABLE `promer` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2016-04-19 12:53:40
