package de.unijena.cheminf.npopensourcecollector.misc;

import org.springframework.stereotype.Service;

import java.util.Arrays;
import java.util.HashSet;


@Service
public class DatabaseTypeChecker {




        private final String[] africa = {"afrodb", "afrocancer", "afromalariadb", "afrotryp", "conmednp", "etm", "mitishamba", "nanpdb", "p-anapl", "sancdb"};
        private final String[] asia = {"him", "hit", "tcmdb_taiwan", "tcmid", "tipdb", "imppat", "inpacdb"};
        private final String[] europe = {"tppt"};
        private final String[] america = {"nubbedb", "uefs", "biofacquim"};

        private final HashSet<String> continentAfrica = new HashSet<String>(Arrays.asList(africa));
        private final HashSet<String> continentAsia = new HashSet<String>(Arrays.asList(asia));
        private final HashSet<String> continentEurope = new HashSet<String>(Arrays.asList(europe));
        private final HashSet<String> continentAmerica = new HashSet<String>(Arrays.asList(america));



        private final String[] plants = {"uefs","tppt","tmdb","tipdb","tcmid", "tcmdb_taiwan","spektraris","sancdb",
                "respect","p-anapl", "npact","nanpdb","mitishamba","inpacdb","imppat", "hitdb","himdb","etmdb","conmednp",
                "afrotryp", "afromalariadb","afrocancer","afrodb", "cmaup"};
        private final String[] bacteria = {"streptomedb"};
        private final String[] fungi = {"lichendatabase"};
        private final String[] animals = {};
        private final String[] marine = {"swmd", "mnp"};
        private final String[] mixed = {"nubbedb","npcare","npatlas", "np_atlas_2019_12","npass","analyticon_all_np", "biofacquim"};

        private final HashSet<String> taxPlants = new HashSet<String>(Arrays.asList(plants));
        private final HashSet<String> taxBacteria = new HashSet<String>(Arrays.asList(bacteria));
        private final HashSet<String> taxFungi = new HashSet<String>(Arrays.asList(fungi));
        private final HashSet<String> taxAnimals = new HashSet<String>(Arrays.asList(animals));
        private final HashSet<String> taxMarine = new HashSet<String>(Arrays.asList(marine));
        private final HashSet<String> taxMixed = new HashSet<String>(Arrays.asList(mixed));


        public  String checkContinent(String sourceDB){

            if(continentAfrica.contains(sourceDB)){
                return "africa";
            }
            else if(continentAsia.contains(sourceDB)){
                return "asia";
            }
            else if(continentEurope.contains(sourceDB)){
                return "europe";
            }
            else if(continentAmerica.contains(sourceDB)){
                return "america";
            }
            else {
                return "nogeo";
            }
        }


        public  String checkKingdom(String sourceDB){
            if(taxPlants.contains(sourceDB)){
                return "plants";
            }
            else if(taxBacteria.contains(sourceDB)){
                return "bacteria";
            }
            else if(taxAnimals.contains(sourceDB)){
                return "animals";
            }
            else if(taxFungi.contains(sourceDB)){
                return "fungi";
            }
            else if(taxMarine.contains(sourceDB)){
                return "marine";
            }
            else if(taxMixed.contains(sourceDB)){
                return "mixed";
            }
            else{
                return "notax";
            }

        }


}
