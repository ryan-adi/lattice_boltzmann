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
                for subtag in subtags:
                    xml_subtag = xml_tag.find(subtag)
                    bc_type = xml_subtag.get("type")
                    if xml_subtag.text==None:
                        self.ctrl_params[tag][subtag] = [bc_type, []]
                    else:
                        attr_data = np.fromstring(xml_subtag.text, dtype=float, sep=" ")
                        self.ctrl_params[tag][subtag] = [bc_type, attr_data]
            else:
                subtags = [child.tag for child in xml_tag]
                for subtag in subtags:
                    xml_subtag = xml_tag.find(subtag)
                    datatype = xml_subtag.get("type")
                    if datatype=="np.array":
                        number_type = xml_subtag.get("number")
                        self.ctrl_params[tag][subtag] = np.fromstring(xml_subtag.text, dtype=number_type, sep=" ")
                    else:    
                        self.ctrl_params[tag][subtag] = eval(datatype)(xml_subtag.text)

        return self.ctrl_params
    
    def get(self):
        return self.ctrl_params
