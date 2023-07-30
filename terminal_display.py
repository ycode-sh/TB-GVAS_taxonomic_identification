#!/bin/python
import pyfiglet


# Make a display

def display_logo():
    logo_text = "dfgmrtc-labs"
    font = "slant"

    # Generate the ASCII art logo
    ascii_art = pyfiglet.figlet_format(logo_text, font=font)

    # Print the logo to the console
    print(ascii_art)

def main():
    #print("Welcome to My Command Line App!")
    display_logo()
    #print("This is a simple example of displaying a logo on the command line.")

if __name__ == "__main__":
    main()
