import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class WebLoader {


    static public void main(String[] args) {
        List<String> urls = new ArrayList<>();
        readUrls(urls, "links.txt");
        WebFrame webFrame = new WebFrame(urls);
    }

    static private void readUrls(List<String> urls, String fileName){
        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(fileName));
            String line;
            while (true) {
                line = reader.readLine();
                if(line == null) break;
                urls.add(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
