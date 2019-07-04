package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.Fragment;
import de.unijena.cheminf.npopensourcecollector.mongocollections.FragmentRepository;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;

@Service
public class FragmentReaderService {

    private LineNumberReader fileReader;

    @Autowired
    public FragmentRepository fragmentRepository;

    public void doWork(Integer withSugar, String file){
        String line;

        try {

            fileReader = new LineNumberReader(new InputStreamReader(new FileInputStream(file)));
            System.out.println("Fragment reader and insertion in MongoDB " + file);


            while ((line = fileReader.readLine()) != null ) {

                String [] tabline = line.split("\\s+"); //0 = signature, 1=score
                Fragment fragment = new Fragment();

                fragment.setHeight(2);
                fragment.setSignature(tabline[0]);
                fragment.setScorenp(Double.parseDouble(tabline[1]));
                fragment.setWith_sugar(withSugar);


                fragmentRepository.save(fragment);
            }
        }catch (IOException ex) {
            System.out.println("Oops ! File not found. Please check if the -in file or -out directory is correct");
            ex.printStackTrace();
            System.exit(0);
        }



    }
}
