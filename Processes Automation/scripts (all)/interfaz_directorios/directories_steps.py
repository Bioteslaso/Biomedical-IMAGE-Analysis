import os
from termcolor import colored
from pathlib import Path



#Function that navigates through files and directories
#x,y are folder names
def files_folders(x, y):

    #Create files and folders empty lists to work with files and folders located in the desired directory
    files = []
    folders = []

    #Join the input folders according to each operating system
    path_element = os.path.join(x, y)

    #List documents inside the path
    documents_inside = os.listdir(path_element)


    #Distinguish between files and folders and append to the corresponding empty list "files" or "folders"
    for element in documents_inside:

        path_file = os.path.join(path_element, element)

        if os.path.isfile(path_file) == True:
            files.append(element)

        elif os.path.isdir(path_file) == True:
            folders.append(element)

    print('\n')

    print('These are the documents inside', colored(path_element, 'green'), 'folder, organized by files and folders:\n ', )


    if len(files) != 0:
        print('Files in this folder are: ')
        for count, element in enumerate(files):
            print(count + 1,'.', element)

    print('\n')


    if len(folders) != 0:
        print('Folders in this folder are: ')
        for count, element in enumerate(folders):
            print(count + 1,'.', element)



    #Boolean to keep asking for the option if the input raises an error

    guessing_document = True

    while guessing_document:

        print('\n')
        print('What do you want to choose?: \n ')
        print('(1) A file')
        print('(2) A folder')
        print('(3) Go back')
        print('')

        try:
            option = input('Choose an option: ')

            if option == '1':

                #Boolean to keep asking for the option if the option introduced is wrong
                guessing_file = True

                while guessing_file:

                    print('')
                    choose_file = input('Choose the number of the file you want to use: ')

                    #Simple way to check if the input is a number
                    numbers = '1234567890'

                    check_numbers = True

                    #If some letter in the input string is not a number, the check_numbes boolean turns False
                    for element in choose_file:
                        if element not in numbers:
                            check_numbers = False

                    #If all the letters in the input string are numbers
                    if check_numbers == True:

                        #Simple way to check if the input is a inside the range of files elements
                        length_files = len(files)
                        length_files_list = list(range(length_files))
                        length_files_list = [str(x) for x in length_files_list]

                        #We have to substract 1 since elements in files list start at 0 index
                        number_file = int(choose_file) - 1

                        if str(number_file) in length_files_list:

                            #Call the chosen file
                            chosen_file = files[number_file]

                            print('')
                            print('You have chosen document: ', colored(chosen_file, 'yellow'))

                            print('Is this correct?')
                            print('(1) Yes')
                            print('(2) No')
                            print('(3) Go back to the current folder menu')

                            #Boolean to keep asking for an input if the introduced one is wrong
                            guessing_option_1 = True

                            while guessing_option_1:

                                option_correct = input("Choose an option: ")

                                if option_correct == '1':
                                    #Path of the file you have chosen
                                    file_path = os.path.join(path_element, chosen_file)
                                    print(colored('The file path is: ', 'green'), colored(file_path, 'yellow'))

                                    #Now we need to terminate all the guessings
                                    #Otherwise, we would get trapped into some of the loops
                                    guessing_option_1 = False
                                    guessing_file = False
                                    guessing_document = False

                                    #End the function
                                    return None

                                elif option_correct == '2':
                                    #If option 2 is requested, the number of the file would be wrong
                                    #Hence we would need to terminate this guessing and go back to the previous one
                                    guessing_option_1 = False

                                elif option_correct == '3':
                                    #If option 3 is selected, the user wants to go back to the previous step, this is, the current folder menu
                                    #Hence, we would need to terminate both this guessing and the previous one
                                    guessing_option_1 = False
                                    guessing_file = False

                                else:
                                    #If neither option 1 or 2 are chosen, then the program asks for a new input
                                    print("You must introduce a correct entry. Try again")

                        else:
                            #If the number file is not in the range of the files list, then the program asks for a new input
                            print("Your entry is not correct. Try again")
                    else:
                        #If the input is not a number or is not well written, then the program asks for a new input
                        print("Your entry is not correct. Try again")

            elif option == '2':

                #To choose the folders option, there must exist folders
                if len(folders) == 0:
                    print(colored("There are no folders in this folder", 'red'))
                    #If there are no folders and you select the option folders, you are redirected to the current folder menu
                else:

                    #Boolean to keep asking for a folder name in case the input is wrong
                    guessing_folder = True

                    while guessing_folder:
                        print('')
                        choose_folder = input('Choose the number of the folder you want to use: ')

                        #Simple way to check if the input is a number
                        numbers = '1234567890'

                        check_numbers = True

                        for element in choose_folder:
                            if element not in numbers:
                                check_numbers = False

                        #If the input introduced is a number
                        if check_numbers == True:

                            length_folders = len(folders)
                            length_folders_list = list(range(length_folders))
                            length_folders_list = [str(x) for x in length_folders_list]

                            number_folder = int(choose_folder) - 1

                            if str(number_folder) in length_folders_list:

                                chosen_folder = folders[number_folder]

                                print('')
                                print('You have chosen folder: ', colored(chosen_folder, 'yellow'))

                                print('Is this correct?')
                                print('(1) Yes')
                                print('(2) No')
                                print('(3) Go back to the current folder menu')

                                guessing_option = True

                                while guessing_option:

                                    option_correct = input("Choose an option: ")

                                    if option_correct == '1':

                                        #Boolean to keep asking for a response in case the input is not correct
                                        guessing_stay = True

                                        #Maybe we want to select a whole folder path or maybe we want to know what is inside it
                                        while guessing_stay:
                                            print('')
                                            print('Now you can choose the correct folder or keep navigating through it:')
                                            print('(1) Choose current folder')
                                            print('(2) Continue navigating through it')

                                            stay_option = input('Choose an option: ')

                                            if stay_option == '1':
                                                #This is the path of the desired folder, it will be printed
                                                path_folder = os.path.join(path_element, chosen_folder)
                                                print(colored('This folder path is: ', 'green'),colored(path_folder, 'yellow'))

                                                #We must terminate all the guessings, since the user decided to choose this path
                                                guessing_option = False
                                                guessing_stay = False
                                                guessing_folder = False
                                                guessing_document = False

                                                #Terminate the function
                                                return

                                            elif stay_option == '2':

                                                #If the user wants to navigate through the folder, all current guessings must be terminated
                                                #The function is know called again but with inputs: current path, folder to which we want to go
                                                guessing_option = False
                                                guessing_stay = False
                                                guessing_folder = False
                                                guessing_document = False
                                                files_folders(path_element, chosen_folder)

                                    elif option_correct == '2':
                                        #If the input is not correct, the program terminates this guessing and asks for a number again
                                        guessing_option = False

                                    elif option_correct == '3':
                                        #If the user wants to go to the current directory menu, guessings must be terminated
                                        guessing_option = False
                                        guessing_folder = False

                                    else:
                                        #If the entry for the option is not correct, the programs asks again for a new one
                                        print('Your entry is not correct. Try again')
                            else:
                                #If the input for the folder number is not correct, the program asks for a new one
                                print("Your entry is not correct. Try again")

                        else:
                            #If the input is not a number, the program asks for a new one
                            print("Your entry is not correct. Try again")

            elif option == '3':
                #If option 3 is selected, the user wants to go to the pevious directory
                one = os.path.abspath(path_element)
                two = os.path.abspath(str(Path.home()))

                if one == two:
                    #Terminates guessing and go to the home directory
                    guessing_document = False
                    files_folders(Path.home(), '')

                else:
                    #Terminates guessing and go to the previous directory
                    guessing_document = False
                    #Get the absolute path
                    x = str(os.path.abspath(path_element))
                    #Split the path in parent path and last directory of the path and select the first element of the list, which is the path
                    x = x.rsplit(os.path.sep, 1)[0]
                    #Call the function with x x = previous path
                    files_folders(x, '')

            else:
                print("Your entry is not correct. Try again")

        except Exception:
            print('Ups')

def directories():

    guessing_directory = True

    print('')
    print('============= DIRECTORY MENU =============')
    print('')
    print('Welcome to the directory menu. Here you can navigate through the directories and choose the target one.\n',)
    print('Bienvenido al menú de directorios. Aquí puedes navegar a través de ellos y elegir aquel sobre el que quieres trabajar.\n')

    #while guessing_directory:

    #Function to get the current directory
    current_directory = os.getcwd()
    print('Your current directory is: ',colored(current_directory, 'green'), '\n')

    print('Now navigate through the directories from your user and choose the path on which you want to work')

    #Get the home path

    user_path = str(Path.home())

    system_separator = os.path.sep

    user_path_list = user_path.split(system_separator)

    #Call the function to display the documents inside the user's folder
    files_folders(user_path, '')



directories()
