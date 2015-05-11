#include "ResultsBuilder.h"

ResultsBuilder::ResultsBuilder (int num_sizes, int num_variations,
                                vector<SizePerm> file_names)
{
    this->num_sizes = num_sizes;
    this->num_variations = num_variations;
    this->file_names = file_names;
}

/* Builds csv tables and write in out_file_name
 *
 */
void ResultsBuilder::buildCSV(vector<Table> tables, string out_file_name)
{
    size_t k;
    ofstream fout;

    fout.open(out_file_name.c_str());

    if ( !fout ) {
        cout << strerror(errno) << '\n'; // displays "Permission denied"
        exit(0);
    }

    for (k = 0; k < tables.size(); k++) {
        Table::iterator it;

        string header = "Solver;";

        OutputLine::iterator it2;

        it = tables[k].begin();
        for (it2 = tables[k][it->F].begin(); it2 != tables[k][it->F].end(); it2++) {
            header+=it2->F+";";
        }

        fout<<header<<endl;

        for (it = tables[k].begin(); it != tables[k].end(); it++) {
            OutputLine::iterator it2;
            string line = it->F+";";
            for (it2 = tables[k][it->F].begin(); it2 != tables[k][it->F].end(); it2++) {
                line += Util::toString((it2->S))+";";
            }
            fout<<line<<endl;

        }
        fout<<endl;
    }

    fout.close();
}

/* Reads the next line in result log file fin. Stores the following:
 * @cost: sorting cost
 * @rev: numver of reversals
 * @seq: sequence of reversals, using the format [[a,b],[c,d]]
 */
void ResultsBuilder::read_result_line(ifstream &fin, ResultLine &result_line)
{
    char buff[MAXBUFF];
    string line;
    getline(fin, line);

    //TODO: use Util::split here
    sscanf(line.c_str(),"%lld %lld %s", &result_line.cost, &result_line.rev,
           buff);
    result_line.seq = string(buff);
}

void ResultsBuilder::write_result_line(ofstream &fout, ResultLine &result_line)
{
    fout<<result_line.cost<<" "<<result_line.rev<<" "<<result_line.seq<<endl;
}

/* Check if result a is lower (e.g better) than b */
bool isLower(ResultLine a, ResultLine b)
{
    if (a.cost < b.cost)
        return true;
    else if (a.cost == b.cost)
        return a.rev < b.rev;
    return false;
}

/* Build a table using same size results of algorithm variations:
 *  @var_fnames: file names  of variations.
 *  @out_fname: file name to store min result
 * */
Table ResultsBuilder::computeTableBySize(vector<Variation> &variations,
        string out_fname)
{
    SequentialDataAnalyser sda;
    Table table;
    int i;
    int num_var = variations.size();
    bool can_read = true;

    /* Stores permutation result of a algorithm variation */
    SourceResult sourceResult;

    /* Input files */
    vector<ifstream *> input(num_var);

    ResultLine result_line;
    /* Min result */
    ResultLine min_result;

    /* Output file */
    ofstream fout(out_fname.c_str());


    /* Open the files */
    for (i = 0; i < num_var; i++) {
        cout<<"Open "<<variations[i].file_name<<endl;
        input[i] = new ifstream(variations[i].file_name.c_str());
    }

    /* Read permutation results */
    while (can_read) { /* For each permutation */
        InstanceResult instanceResult;
        /* Initialize minimum result for current permutation */
        min_result.cost = INF;
        min_result.rev = INF;

        /* Read results of each algorithm variation */
        for (i = 0; i < num_var; i++) {
            read_result_line(*input[i], result_line);

            /* Check for end of input */
            if (input[i]->eof()) {
                cout<<"Can't read"<<endl;
                can_read = false;
                break;
            }

            /* Fill the results fields */
            sourceResult.name = variations[i].name;
            sourceResult.field["cost"] = result_line.cost;
            sourceResult.field["rev"] = result_line.rev;
            instanceResult.push_back(sourceResult);

            /* Compute minimum result */
            if (isLower(result_line, min_result)) {
                min_result = result_line;
            }
        }

        if (can_read) {
            /* Write minimun result */
            write_result_line(fout, min_result);

            /* Fill the results fields */
            sourceResult.name = "min_result";
            sourceResult.field["cost"] = min_result.cost;
            sourceResult.field["rev"] = min_result.rev;
            instanceResult.push_back(sourceResult);
            /* Add instance (permutation) results to analyser */
            sda.addIns(instanceResult);
        }
    }

    /* Close the opened files */
    for (i = 0; i < num_var; i++) {
        input[i]->close();
    }

    /* Make table */
    table = sda.makeTable();

    return table;
}
