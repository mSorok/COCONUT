package de.unijena.cheminf.npopensourcecollector.misc;

import org.springframework.stereotype.Service;

import java.util.Arrays;
import java.util.HashSet;


@Service
public class DatabaseTypeChecker {




        private final String[] africa = {"afrodb", "afrocancer", "afromalariadb", "afrotryp", "conmednp", "etm", "mitishamba", "nanpdb", "p-anapl", "sancdb"};
        private final String[] china = {"him", "hit", "tcmdb_taiwan", "tcmid", "tipdb"};
        private final String[] india = {"imppat", "inpacdb"};
        private final String[] europe = {"tppt"};
        private final String[] america = {"nubbedb", "uefs", "biofacquim"};

        private final HashSet<String> continentAfrica = new HashSet<String>(Arrays.asList(africa));
        private final HashSet<String> continentIndia = new HashSet<String>(Arrays.asList(india));
        private final HashSet<String> continentChina = new HashSet<String>(Arrays.asList(china));
        private final HashSet<String> continentEurope = new HashSet<String>(Arrays.asList(europe));
        private final HashSet<String> continentAmerica = new HashSet<String>(Arrays.asList(america));



        private final String[] plants = {"uefs","tppt","tmdb","tipdb","tcmid", "tcmdb_taiwan","spektraris","sancdb",
                "respect","p-anapl", "npact","nanpdb","mitishamba","inpacdb","imppat", "hit","him","etm","conmednp",
                "afrotryp", "afromalariadb","afrocancer","afrodb"};
        private final String[] bacteria = {"streptomedb"};
        private final String[] fungi = {"lichendatabase"};
        private final String[] animals = {};
        private final String[] marine = {};
        private final String[] mixed = {"nubbedb","npcare","npatlas","npass","analyticon_all_np", "biofacquim"};

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
            else if(continentChina.contains(sourceDB)){
                return "china";
            }
            else if(continentIndia.contains(sourceDB)){
                return "india";
            }
            else if(continentEurope.contains(sourceDB)){
                return "europe";
            }
            else if(continentAmerica.contains(sourceDB)){
                return "southamerica";
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
