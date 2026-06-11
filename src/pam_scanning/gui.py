"""Tkinter GUI front-end for PAM-scanning.

Collects run parameters through a simple form and hands them to
:func:`pam_scanning.chimeras.pamscan`. Launch with ``pam-scan-gui``.
"""

import tkinter as tk
from tkinter import filedialog, Toplevel, scrolledtext, messagebox
from os import sep, getcwd


def main():

    # Function to show a pop-up window with a text box for information
    def show_info(info_message):
        info_window = Toplevel(root)
        info_window.title("Information")
        info_window.geometry("400x200")
    
        # Create a scrolled text box widget inside the pop-up window
        info_text = scrolledtext.ScrolledText(info_window, wrap=tk.WORD, width=50, height=8)
        info_text.pack(padx=10, pady=10)
    
        # Insert the description text into the text box
        info_text.insert(tk.END, info_message)
        info_text.config(state=tk.DISABLED)  # Make the text box read-only

    # Function to open file dialog and display the selected file path
    def open_file(file_label):
        file_path = filedialog.askopenfilename(initialdir=getcwd())  # Default directory set to root
        if file_path:
            file_label.config(text=f"Selected: {file_path}")
            print(f"File path: {file_path}")
        return file_path

    # Function to open directory dialog and display the selected directory path
    def open_directory(directory_label):
        directory_path = filedialog.askdirectory(initialdir=getcwd())  # Set the default directory to the root directory
        if directory_path:
            directory_label.config(text=f"Save Directory: {directory_path}")
            print(f"Directory path: {directory_path}")
        return directory_path

    # Function to validate integer input
    def validate_integer_input(new_value):
        if new_value.isdigit() or new_value == "":
            return True
        return False

    # Function to retrieve the value of a string entry
    def get_string_value(entry):
        return entry.get()

    # Function to retrieve the value of an integer entry
    def get_integer_value(entry):
        try:
            return int(entry.get())
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter a valid integer.")
            return None

    # Create the main window
    root = tk.Tk()
    root.title("File Open, String, and Integer Variables")
    root.geometry("1200x900")

    # Define fixed column widths
    label_width = 50
    button_width = 50
    info_button_width = 10

    # File Open Dialogs (Rows 1-3)
    file_open_dialogs = ["Open File: ORF", "Open File: ORF+", "Open File: Genome Sequence", "Open File: Codon Table", "Open File: Codon Selection File"]
    file_info = [
                    "Open the open reading frame (ORF) fasta formatted file that contains the gene of interest. " +\
                    "The ORF sequence should begin with the ATG start codon and end with a stop codon.",

                    "Open the open reading frame (ORF) fasta formatted file that contains the gene of interest flanked by 1000 bp of genome homology.",

                    "The host genome sequence in fasta format. For example, if the PAM scan is to be done in yeast, this is the full yeast genome sequence.",

                    "The codon table for the host genome sequence. For formatting, refer to the yeast codon table file in the docs: yeast_64_1_1_all_nuclear.cusp.txt.",

                    "A codon selection file indicating specific chimera insertion points for the PAM scan. " +\
                    "Providing this file overrides the codon sampling frequency parameter below."
    ]
    file_labels = []
    for i in range(5):
        # file_label = tk.Label(root, text="No file selected", width=label_width, anchor='w', wraplength=450)
        file_label = tk.Label(root, text="No file selected", width=label_width, anchor='w', wraplength=450)
        file_label.grid(row=i, column=0, padx=10, pady=10)
        file_labels.append(file_label)

        file_button = tk.Button(root, text=file_open_dialogs[i], width=button_width, command=lambda l=file_label: open_file(l))
        file_button.grid(row=i, column=1, padx=10, pady=10)

        info_button = tk.Button(root, text="Info", width=info_button_width, command=lambda i=i: show_info(file_info[i]))
        info_button.grid(row=i, column=2, padx=10, pady=10)

    # String Variables (Rows 6-10)
    string_labels = ["Gene name", "Local BLAST database", "Guide primer: forward suffix", "Insert primer: forward suffix", "Insert primer: reverse suffix"]
    string_entries_default = ["MFG", "yeast", "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAG", "GAAGATGTTGTCTGTTGCTCTATGTCATAT", "CTTCTACAACAGACAACGAGATACAGTATA"]
    string_info = [

                    "Provide a gene name label for the calculation.",

                    "The name of your local BLAST database. To create a local BLAST database, you must first install BLAST+: " +\
                    "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. " +\
                    "There are several good BLAST+ tutorials available online.",

                    "This forward primer sequence amplifies the yeast CRISPR plasmid starting at the base after the unique targeting sequence (UnTS). " +\
                    "This sequence will be prepended with the 20 bp guide sequence needed for each PAM scan chimera site. " +\
                    "The composite primer is then used to build a complete guide plasmid using our around-the-world protocol. ",

                    "The sequence for the forward top strand, 5' to 3', insertion primer to be concatenated downstream of 5' to 3' PAM scan chimera site genome homology. " +\
                    "This forward primer, together with the reverse insertion primer below, is used to amplify the chimeric insert for each specific PAM scan site (i.e., codon)",

                    "The sequence for the reverse bottom strand, 3' to 5', insertion primer to be concatenated upstream of the PAM scan chimera site genome homology. " +\
                    "This reverse primer, together with the forward insertion primer above, is used to amplify the chimeric insert for each specific PAM scan site (i.e., codon)."

    ]
    string_entries = []
    for i in range(5):
        string_label = tk.Label(root, text=string_labels[i], width=label_width, anchor='w')
        string_label.grid(row=i+5, column=0, padx=10, pady=5)

        string_var = tk.StringVar(value=string_entries_default[i])
        string_entry = tk.Entry(root, text=string_var, width=button_width)
        string_entry.grid(row=i+5, column=1, padx=10, pady=5)
        string_entries.append(string_entry)

        info_button = tk.Button(root, text="Info", width=info_button_width, command=lambda i=i: show_info(string_info[i]))
        info_button.grid(row=i+5, column=2, padx=10, pady=5)

    # Validation command for integer inputs
    vcmd = (root.register(validate_integer_input), '%P')

    # Integer Variables (Rows 10-13)
    int_labels = ["Primer length (bp)", "Max gap between two PAMs (bp)", "Codon sampling frequency", "Max PAM inclusions", "Max PAM inclusion length (bp)"]
    integer_entries_default = [100, 60, 1, 5, 15]
    int_info = [
                    "Set the length of DNA primers that will be used to amplify the chimera insert.",

                    "Sequential PAM sites separated by more than this distance are expected to create PAM scanning gaps " +\
                    "anticipated to have lower editing efficiency. 60 is the default value that reflects the empirical 30 bp rule: " +\
                    "chimeric insertion at codons >30 bp from the CRISPR cut sight is less efficient.",

                    "With the exception of PAM scanning gaps, this parameter sets the codons at which the chimera sequence inserted. " +\
                    "A frequency of 1, would be exhaustively after every codon. A frequence of 2 corresponds to every other codon, and so on. " +\
                    "This parameter is ignored if a codon selection file is provided.",

                    "Guide sequences are defined as 20 bp in length. However, slightly shorter guide sequences may still enable cutting. " +\
                    "We call these PAM inclusions and we track them across all guide solutions. 5 is the default value because it generates ample guide designs " +\
                    "with a minimal number of PAM inclusions. The idea is to design a guide with minimal PAM inclusions " +\
                    "and sequence gaps between the Cas9 cut sight and codon insertion point as enforced by a 30 bp rule. " +\
                    "Providing this file overrides the codon sampling frequency parameter below.",

                    "Guide sequences are defined as 20 bp in length. However, slightly shorter guide sequences may still enable cutting. " +\
                    "We call these PAM inclusions and we track them across all guide solutions. This parameter defines the minimum length of a PAM inclusion. " +\
                    "Using the default value of 15 as an example, all candidate PAM inclusions having length > 15 are counted as PAM inclusions. "

    ]
    integer_entries = []
    for i in range(5):
        int_label = tk.Label(root, text=int_labels[i], width=label_width, anchor='w')
        int_label.grid(row=i+11, column=0, padx=10, pady=5)

        int_var = tk.IntVar(value=integer_entries_default[i])
        int_entry = tk.Entry(root, textvariable=int_var, width=button_width, validate="key", validatecommand=vcmd)
        int_entry.grid(row=i+11, column=1, padx=10, pady=5)
        integer_entries.append(int_entry)

        info_button = tk.Button(root, text="Info", width=info_button_width, command=lambda i=i: show_info(int_info[i]))
        info_button.grid(row=i+11, column=2, padx=10, pady=5)

    # Save Directory Row (Row 14)
    save_label = tk.Label(root, text="No directory selected", width=label_width, anchor='w', wraplength=450)
    save_label.grid(row=16, column=0, padx=10, pady=10)

    save_button = tk.Button(root, text="Select Save Directory", width=button_width, command=lambda: open_directory(save_label))
    save_button.grid(row=16, column=1, padx=10, pady=10)

    info_button = tk.Button(root, text="info", width=info_button_width, command=lambda: show_info("Select a save directory for the PAM scan output files"))
    info_button.grid(row=16, column=2, padx=10, pady=10)

    # Function to print out values
    def print_values():

        file_paths = {  
                        'orf_file_path':"", 
                        'orf_plus_buffer_file_path':"", 
                        'local_genome_file_path':"", 
                        'codon_table_file_path':"",
                        'codon_selection_file_path':""
                        
                        }

        # print("File paths selected:")
        for i in range(5):
            text = file_labels[i].cget('text')
            if text != "No file selected":
                text = text.split("Selected: ")[1]
            file_paths[list(file_paths)[i]] = text
            # print(file_labels[i].cget('text'))
            # print(f"File {i+1}: {file_labels[i].cget('text')}")
        # orf_file_path = file_paths['orf_file_path']
        # orf_plus_buffer_file_path = file_paths['orf_plus_buffer_file_path']
        # local_genome_file_path = file_paths['local_genome_file_path']
        # codon_selection_file_path = file_paths['codon_selection_file_path']

        string_variables = {  
                            'geneName':"", 
                            'localBlastDb':"", 
                            'guidePrimerForwardSuffix':"", 
                            'insertPrimerForwardSuffix':"",
                            'insertPrimerReverseSuffix':""
                        
                            }
    
        # print("\nString variable values:")
        for i, entry in enumerate(string_entries):
            string_variables[list(string_variables)[i]] = get_string_value(entry)
            # print(f"String Variable {i+1}: {get_string_value(entry)}")
        # geneName = string_variables['geneName']
        # localBlastDb = string_variables['localBlastDb']
        # guidePrimerForwardSuffix = string_variables['guidePrimerForwardSuffix']
        # insertPrimerForwardSuffix = string_variables['insertPrimerForwardSuffix']
        # insertPrimerReverseSuffix = string_variables['insertPrimerReverseSuffix']

        int_variables = {  
                            'primerLength':0, 
                            'maxPamCutGap':0, 
                            'codonsSamplingGap':0, 
                            'pamInclusionThreshold':0,
                            'pamInclusionSequenceThreshold':0,
                            'codonsSamplingGap':0

                            }

        # print("\nInteger variable values:")
        for i, entry in enumerate(integer_entries):
            int_variables[list(int_variables)[i]] = get_integer_value(entry)
            # print(f"Integer Variable {i+1}: {get_integer_value(entry)}")
        # primerLength = int_variables['primerLength']
        # maxPamCutGap = int_variables['maxPamCutGap']
        # codonsSamplingGap = int_variables['codonsSamplingGap']
        # pamInclusionThreshold = int_variables['pamInclusionThreshold']
        # pamInclusionSequenceThreshold = int_variables['pamInclusionSequenceThreshold'] 
        # codonsSamplingGap = int_variables['codonsSamplingGap']

        text = save_label.cget('text')
        if text != "No file selected":
            text = text.split("Save Directory: ")[1]
        save_locations = {"outputPath":text}
        # print(f"\nSave location: {save_label.cget('text')}")

        from pam_scanning.chimeras import pamscan
        kwargs = file_paths
        kwargs.update(string_variables)
        kwargs.update(int_variables)
        kwargs.update(save_locations)
        pamscan(**kwargs)

    # Button to print out all values
    submit_button = tk.Button(root, text="submit", command=print_values)
    submit_button.grid(row=17, column=1, padx=10, pady=20)

    # Start the GUI event loop
    root.mainloop()

if __name__ == "__main__":
    main()
