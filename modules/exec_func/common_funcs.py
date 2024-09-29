from textual.widgets import Input


def get_input(self, input_widget_id):
    return str(self.query_one(f"#{input_widget_id}", Input).value).strip()
