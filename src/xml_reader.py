from common_modules import ET, np

class XMLReader(): 
    def __init__(self):
        self.ctrl_params = {}

    def read(self, xml_file):
        tree = ET.parse(xml_file)
        xml_root = tree.getroot()

        self.ctrl_params["caseName"] = xml_root.get("caseName")
        tags = [child.tag for child in xml_root]

        for tag in tags:
            self.ctrl_params[tag] = {}
            xml_tag = xml_root.find(tag)

            if tag=="BoundaryConditions":
                subtags = [child.tag for child in xml_tag]
                if "corner" in subtags:
                    self.ctrl_params[tag]["corner"] = int(xml_tag.find("corner").text)
                else:
                    self.ctrl_params[tag]["corner"] = 0

                bc_loc = ["left", "right", "top", "bottom"]
                for subtag in bc_loc:
                    xml_subtag = xml_tag.find(subtag)
                    bc_type = xml_subtag.get("type")
                    if xml_subtag.text==None:
                        self.ctrl_params[tag][subtag] = [bc_type, []]
                    else:
                        attr_data = np.fromstring(xml_subtag.text, dtype=float, sep=" ")
                        self.ctrl_params[tag][subtag] = [bc_type, attr_data]
            elif tag=="Particles":
                tag_iter = xml_tag.iter()
                next(tag_iter)
                for elem in tag_iter:
                    xml_subtag = elem
                    subtag = elem.tag

                    datatype = xml_subtag.get("type")
                    number_type = xml_subtag.get("number")
                    if subtag not in self.ctrl_params[tag]:
                        first = np.fromstring(xml_subtag.text, dtype=number_type, sep=" ")
                        self.ctrl_params[tag][subtag] = np.stack([first], axis=0)
                    else:
                        first = self.ctrl_params[tag][subtag]
                        second = np.fromstring(xml_subtag.text, dtype=number_type, sep=" ")
                        second = np.stack([second], axis=0)
                        self.ctrl_params[tag][subtag] = np.concatenate((first, second), axis=0)
            else:
                tag_iter = xml_tag.iter()
                next(tag_iter)
                for elem in tag_iter:
                    xml_subtag = elem
                    subtag = elem.tag

                    datatype = xml_subtag.get("type")
                    if datatype=="np.array":
                        number_type = xml_subtag.get("number")
                        self.ctrl_params[tag][subtag] = np.fromstring(xml_subtag.text, dtype=number_type, sep=" ")
                    else:    
                        self.ctrl_params[tag][subtag] = eval(datatype)(xml_subtag.text)

        return self.ctrl_params
    
    def get(self):
        return self.ctrl_params
